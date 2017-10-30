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

Spectrum MicroBsdf::Sample_f(const Vector3f &woWorld, Vector3f *wiWorld,  Sampler &sampler,
                      Float *pdf, BxDFType type, BxDFType *sampledType) const
{
    Spectrum f (1);
    Vector3f wi, wo = WorldToLocal(woWorld);

    ProfilePhase pp(Prof::BSDFSampling);
    // Choose which _BxDF_ to sample
    int matchingComps = NumComponents(type);
    BxDF *bxdf = bxdfs[0];
    MicrofacetReflection *mr = (MicrofacetReflection*) bxdf;
    const MicrofacetDistribution *dist = mr->getDistribution();

    Float phi, theta;
    Vector3f wh = dist->Sample_wh(wo, sampler.Get2D(), &phi, &theta);
    *pdf = dist->D(wh);

    Groove groove(phi, theta);
    int bounce = groove.bounce(wo, sampler.Get1D(), wi);
    if (bounce > 0 && bounce < groove.MAX_BOUNCE) {
        *wiWorld = LocalToWorld(wi);
        if (bounce > 1) printf("bounce count: %d\n", bounce);
        return f= *pdf * bounce;
        //return bounce;
    } else {
        //in shadow?
        *pdf = 0;
        return 0;
    }
}

}//namespace

