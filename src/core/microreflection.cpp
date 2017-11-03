#include "microreflection.h"
#include "microfacet.h"
#include "spectrum.h"
#include "sampler.h"
#include "sampling.h"
#include "interpolation.h"
#include "scene.h"
#include "interaction.h"
#include "stats.h"
#include <stdarg.h>

namespace pbrt{

Spectrum MicroBsdf::f(const Vector3f &woW, const Vector3f &wiW,
                 BxDFType flags) const {
    return 0;
}

Float MicroBsdf::Pdf(const Vector3f &woWorld, const Vector3f &wiWorld,
                BxDFType flags) const {
    return 0;
}

Spectrum MicroBsdf::Sample_f(const Vector3f &woWorld, Vector3f *wiWorld,  Sampler &sampler,
                      Float *pdf, BxDFType type, BxDFType *sampledType) const
{
    Spectrum f (1);
    Vector3f wi, wo = WorldToLocal(woWorld);

    ProfilePhase pp(Prof::BSDFSampling);
    // Choose which _BxDF_ to sample
    //int matchingComps = NumComponents(type);
    BxDF *bxdf = bxdfs[0];
    MicrofacetReflection *mr = (MicrofacetReflection*) bxdf;
    const MicrofacetDistribution *dist = mr->getDistribution();

    Float phi, theta;
    Vector3f wh = dist->Sample_wh(wo, sampler.Get2D(), &phi, &theta);
    Float cosThetaH = AbsDot(wo, wh);

    *pdf = dist->Pdf(wo, wh) / (4 * cosThetaH);
    //*pdf = dist->Pdf(wo, wh);
    //*pdf = 1;


    Groove groove(phi, theta);
    float bounceCount = groove.bounce(wo, sampler.Get1D(), wi);
    if (bounceCount > 0 && bounceCount < groove.MAX_BOUNCE && SameHemisphere(wo, wi)) {
       
        Float cosThetaO = AbsCosTheta(wo), cosThetaI = AbsCosTheta(wi);
        if (bounceCount > 1) printf("bounce count: %f\n", bounceCount);
        *wiWorld = LocalToWorld(wi);
        //return f= dist->D(wh) * bounceCount / (4.0 * cosThetaO * cosThetaI) ;
        return f= dist->D(wh)  / (4.0 * cosThetaO * cosThetaI) ;
        //return f=  bounceCount / (4.0 * cosThetaO * cosThetaI) ;
        //return f = bounceCount;
    } else {
        //in shadow?
        *pdf = 0;
        return 0;
    }
}

}//namespace

