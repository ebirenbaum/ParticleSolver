#include "gasconstraint.h"

GasConstraint::GasConstraint(double density, QList<int> *particles)
    : Constraint(), p0(density)
{
    neighbors = new QList<int>[particles->size()];
    deltas = new glm::dvec2[particles->size()];

    for (int i = 0; i < particles->size(); i++) {
        ps.append(particles->at(i));
    }
}

GasConstraint::~GasConstraint()
{
    delete[] neighbors;
}

void GasConstraint::project(QList<Particle *> *estimates)
{
    // Find neighboring particles and estimate pi for each particle
    lambdas.clear();
    for (int k = 0; k < ps.size(); k++) {
        neighbors[k].clear();
        int i = ps[k];
        Particle *p_i = estimates->at(i);
        double pi = 0., denom = 0.;

        // Find neighbors
        for (int j = 0; j < estimates->size(); j++) {

            // Check if the next particle is actually this particle
            if (j != i) {
                Particle *p_j = estimates->at(j);
                glm::dvec2 r = p_i->ep - p_j->ep;
                double rlen2 = glm::dot(r, r);
                if (rlen2 < H2) {

                    // Found a neighbor! Remember it and add to pi and the gamma denominator
                    neighbors[k].append(j);
                    double incr = poly6(rlen2) / p_j->imass;
                    if (p_j->ph == SOLID) {
                        incr *= S_SOLID;
                    }
                    pi += incr;

                    glm::dvec2 gr = grad(estimates, k, j);
                    denom += glm::dot(gr, gr);
                }

            // If it is, cut to the chase
            } else {
                neighbors[k].append(j);
                pi += poly6(0) / p_i->imass;
            }
        }

        glm::dvec2 gr = grad(estimates, k, i);
        denom += glm::dot(gr, gr);

        // Compute the gamma value
//        cout << i << " estimated " << pi << endl;
        double p_rat = (pi/p0);
        if(p_rat < 1) p_rat = 1;
        double lambda = -(p_rat - 1.) / (denom + RELAXATION);
        lambdas[i] = lambda;
    }

    // Compute actual deltas
    for (int k = 0; k < ps.size(); k++) {
        glm::dvec2 delta = glm::dvec2();
        int i = ps[k];
        Particle *p_i = estimates->at(i);

        for (int x = 0; x < neighbors[k].size(); x++) {
            int j = neighbors[k][x];
            if (i == j) continue;
            Particle *p_j = estimates->at(j);
            glm::dvec2 r = p_i->ep - p_j->ep;
            double rlen = glm::length(r);
            glm::dvec2 sg = spikyGrad(r, rlen);
            double lambdaCorr = -K_P * pow((poly6(rlen) / poly6(DQ_P * H)), E_P);
            delta += (lambdas[i] + lambdas[j] + lambdaCorr) * sg;
        }
        deltas[k] = (delta / p0);
    }

    for (int k = 0; k < ps.size(); k++) {
        int i = ps[k];
        Particle *p_i = estimates->at(i);
        p_i->ep += deltas[k] / (double) neighbors[k].size();
    }
}

void GasConstraint::draw(QList<Particle *> *particles)
{

}

double GasConstraint::poly6(double r2)
{
    if(r2 >= H2) return 0;
    double term2 = (H2 - r2);
    return (315. / (64. * M_PI * H9)) * (term2 * term2 * term2);
//    return (H-r) / (H*H);
}

glm::dvec2 GasConstraint::spikyGrad(const glm::dvec2 &r, double rlen)
{
    if(rlen >= H) return glm::dvec2();
    if(rlen == 0) return glm::dvec2();
    return -glm::normalize(r) * (45. / (M_PI * H6)) * (H - rlen) * (H - rlen);
//    return -r / (H*H*rlen);
}

glm::dvec2 GasConstraint::grad(QList<Particle *> *estimates, int k, int j)
{
    int i = ps[k];
    Particle *p_i = estimates->at(i), *p_j = estimates->at(j);
    glm::dvec2 r = p_i->ep - p_j->ep;
    double rlen = glm::length(r);
    if (p_i != p_j) {
        return -spikyGrad(r, rlen) / (p0);
    }

    glm::dvec2 out = glm::dvec2();
    for (int x = 0; x < neighbors[k].size(); x++) {
        r = p_i->ep - estimates->at(neighbors[k][x])->ep;
        rlen = glm::length(r);
        out += spikyGrad(r, rlen);
    }

    return out / (p0);
}

double GasConstraint::evaluate(QList<Particle *> *estimates)
{
    std::cout << "You shouldn't be calling evaluate on fluids" << std::endl;
    exit(1);
}

glm::dvec2 GasConstraint::gradient(QList<Particle *> *estimates, int respect)
{
    std::cout << "You shouldn't be calling gradient on fluids" << std::endl;
    exit(1);
}

void GasConstraint::updateCounts(int *counts)
{
    std::cout << "You shouldn't be calling update counts on fluids yet" << std::endl;
    exit(1);
//    for(int i = 0; i < numParticles; i++) {
//        Particle *p = particles[i];
//    }
//    counts[idx]++;
}
