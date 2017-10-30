
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_MICROREFLECTION_H
#define PBRT_CORE_MICROREFLECTION_H

// core/reflection.h*
#include "reflection.h"
#include "microfacet.h"

namespace pbrt {

#define SMALLNUMBER 1e-6

class GroovePlane {
public:
    GroovePlane() {};
    GroovePlane(const Vector3f& n, Float xb, Float yb, int xsign): mNormal(n), mXb(xb), mYb(yb), mXsign(xsign)
    {
        mPoint.x = mXb * mXsign;
        mPoint.y = 0;
        mPoint.z = 0;

    };
    Vector3f mNormal;
    Point3f mPoint;
    Float mXb, mYb;
    int mXsign;
    
    Float trace(const Ray& ray, Point3f& ip, bool& degenerate) 
    {
        Float Nx = mNormal[0] * mXb * mXsign;
        Float dDotN = Dot(ray.d, mNormal);
        //Vector3f ro (ray.o.x, ray.o.y, ray.o.z);
        Float oDisP = Dot(mNormal, mPoint - ray.o);
        if (fabs(oDisP) < SMALLNUMBER) {
            degenerate = true;
            return -1;
        }
        if (fabs(dDotN) > SMALLNUMBER){
            Float t = oDisP/dDotN;
            if (t > 0) {
                ip = ray.o + ray.d * t;
                if (ip[1] > -mYb && ip[1] < mYb) {
                    if (mXsign == 1) {
                        if (ip[0] < mXb && ip[0] > 0) {
                            return t;
                        }
                    } else {
                        if (ip[0] > -mXb && ip[0] < 0) {
                            return t;
                        }
                    }
                }
            }
        }
        return -1;
    }
};

class Groove{
public:

    Groove(Float phi, Float thetah) {
        //RADIUS = .01;
        RADIUS = 1;
        MAX_BOUNCE = 5;
        xradius = RADIUS * 10;
        yradius = RADIUS * 1000;
        Transform t1 = RotateY(thetah);
        Transform t2 = RotateY(-thetah);
        Vector3f N1 = t1(Vector3f(0, 0, 1));
        Vector3f N2 = t2(Vector3f(0, 0, 1));

        //Feng to test phi at 0

        //phi = 0;

        grooveXform = RotateZ(phi);
        grooveInvXform = RotateZ(-phi);
        groove_ss = grooveXform(Vector3f(1, 0, 0));
        groove_ts = grooveXform(Vector3f(0, 1, 0));
        groove_ns = grooveXform(Vector3f(0, 0, 1));  

        
        groove_N1 = grooveXform(N1);
        groove_N2 = grooveXform(N2);
        planes[0] = GroovePlane(groove_N2, xradius, yradius, 1);
        planes[1] = GroovePlane(groove_N1, xradius, yradius, -1);

    }
    Vector3f LocalToGroove(const Vector3f &v) const {
        return grooveInvXform(v);
        //return Vector3f(Dot(v, groove_ss), Dot(v, groove_ts), Dot(v, groove_ns));
    }
    Point3f LocalToGroove(const Point3f &p) const {
        return grooveInvXform(p);
        //return Vector3f(Dot(v, groove_ss), Dot(v, groove_ts), Dot(v, groove_ns));
    }


    Vector3f GrooveToLocal(const Vector3f& v) const{
        /*
        return Vector3f(groove_ss.x * v.x + groove_ts.x * v.y + groove_ns.x * v.z,
                            groove_ss.y * v.x + groove_ts.y * v.y + groove_ns.y * v.z,
                                                    groove_ss.z * v.x + groove_ts.z * v.y + groove_ns.z * v.z);
        */
        return grooveXform(v);
    }


    void createRay(const Vector3f& wo, Float u, Ray& ray)
    {
        Vector3f grooveWo = LocalToGroove(wo);
        Point3f ro(u * xradius * 2 - xradius, 0, 0);
        Point3f grooveO = LocalToGroove(ro);
        grooveO = grooveO + grooveWo * (.001); //bias the ray orign back along the ray direction
        ray.o = grooveO;
        ray.d = -grooveWo;
    }

    int otherPlaneIndex(int i) 
    {
        return i == 0 ? 1 : 0;
    }

    int bounce(const Vector3f& wo, const Float u, Vector3f& wi) {

        Ray ray;
        createRay(wo, u, ray);
       
        bool hit = false;
        Point3f hitPt;
        Vector3f refD;
        Ray reflectRay;
        int bounceCount = 0;
        int curHitPlaneIndex = -1;
        bool degenerate = false;

        for (int i = 0; i < 2; i++) {
            if (planes[i].trace(ray, hitPt, degenerate) > 0) {
                hit = true;
                refD = Reflect(-ray.d, planes[i].mNormal);
                reflectRay.o = hitPt;
                reflectRay.d = refD;
                bounceCount++;
                curHitPlaneIndex = i;
                break;
            }
        }
        while (hit && bounceCount < MAX_BOUNCE) {
            int pIdx = otherPlaneIndex(curHitPlaneIndex);
            Vector3f diff = hitPt - planes[pIdx].mPoint;
            Float pside = Dot(diff, planes[pIdx].mNormal);
            Float rayside = Dot(reflectRay.d, planes[pIdx].mNormal); 
            if (pside * rayside > 0) {
                Vector3f groove_wi = reflectRay.d;
                wi = GrooveToLocal(groove_wi);
                hit = false;
                break;
            }
            if (planes[pIdx].trace(reflectRay, hitPt, degenerate) > 0) {
                refD = Reflect(-reflectRay.d, planes[pIdx].mNormal);
                reflectRay.o = hitPt;
                reflectRay.d = refD;
                bounceCount++;
                hit = true;
                curHitPlaneIndex = pIdx;
            } else {
                Vector3f groove_wi = reflectRay.d;
                wi = GrooveToLocal(groove_wi);
                hit = false;
                break;
            }
        }
        return bounceCount;
            
    }

    GroovePlane planes[2];
    Float RADIUS, xradius, yradius; 
    Vector3f groove_ss, groove_ts, groove_ns, groove_N1, groove_N2;
    Transform grooveXform, grooveInvXform;
    int MAX_BOUNCE;

};



class MicroBsdf : public BSDF {
public:
    MicroBsdf(const SurfaceInteraction &si, Float eta = 1): BSDF(si, eta) {}
    virtual Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, Sampler &sampler,
                      Float *pdf, BxDFType type = BSDF_ALL,
                      BxDFType *sampledType = nullptr) const;
    /*
    virtual Spectrum f(const Vector3f &woW, const Vector3f &wiW,
               BxDFType flags = BSDF_ALL) const;
    virtual Float Pdf(const Vector3f &wo, const Vector3f &wi,
              BxDFType flags = BSDF_ALL) const;
    virtual Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType type = BSDF_ALL,
                      BxDFType *sampledType = nullptr) const;
    */

private:
};

}  // namespace pbrt

#endif  // PBRT_CORE_REFLECTION_H
