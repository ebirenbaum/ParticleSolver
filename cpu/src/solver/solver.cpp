#include "solver.h"

Solver::Solver()
{
    m_b = new double[2];
    m_gamma = new double[2];
    m_nCons = -1;

    m_dp = new double[2];
    m_counts = new int[2];
    m_nParts = -1;
}

Solver::~Solver()
{
    delete[] m_b;
    delete[] m_gamma;
    delete[] m_dp;
    delete[] m_counts;
}

int Solver::getCount(int idx)
{
    return m_counts[idx];
}

void Solver::setupM(QList<Particle *> *particles, bool contact)
{
    m_invM.reset(particles->size() * 2, particles->size() * 2);
    for (int i = 0; i < particles->size(); i++) {

        // Down the diagonal
        Particle *p = particles->at(i);
        m_invM.setValue(2 * i, 2 * i, contact ? p->tmass : p->imass);
        m_invM.setValue(2 * i + 1, 2 * i + 1, contact ? p->tmass : p->imass);
    }
}

void Solver::setupSizes(int numParts, QList<Constraint *> *constraints)
{
    bool change = true;
    int numCons = constraints->size();

    // Only update some things if the number of particles changed
    if (m_nParts != numParts) {
        m_nParts = numParts;
        delete[] m_dp;
        delete[] m_counts;
        m_dp = new double[m_nParts * 2];
        m_counts = new int[m_nParts];
        change = true;
    }

    // Update how many constraints affect each particle
    for (int i = 0; i < m_nParts; i++) {
        m_counts[i] = 0;
    }
    for (int i = 0; i < numCons; i++) {
        constraints->at(i)->updateCounts(m_counts);
    }

    if (m_nCons != numCons) {
        m_nCons = numCons;
        delete[] m_b;
        delete[] m_gamma;
        m_b = new double[m_nCons];
        m_gamma = new double[m_nCons];
        change = true;
    }

    if (change) {
        m_JT.reset(m_nParts * 2, m_nCons);
    }
}

void Solver::solveAndUpdate(QList<Particle *> *particles, QList<Constraint *> *constraints, bool stabile)
{
    if (constraints->size() == 0) {
        return;
    }

    // Reset J^T and b
    bool updated = false;

    // Loop!
    for (int i = 0; i < particles->size(); i++) {
        for (int j = 0; j < constraints->size(); j++) {
            Constraint *cons = constraints->at(j);

            // Update b
            if (!updated) {
                m_b[j] = -cons->evaluate(particles);
            }
            glm::vec2 grad_ji = cons->gradient(particles, i);
            m_JT.setValue(2 * i, j, grad_ji.x);
            m_JT.setValue(2 * i + 1, j, grad_ji.y);
        }
        updated = true;
    }

    SparseMatrix temp = m_invM * m_JT;
    m_A = m_JT.getTranspose() * temp;
//    m_JT.printMatrix(4,false);
    m_eq.setA(&m_A);
//    cout << endl;
//    for (int i = 0; i < particles->size(); i++) {
//        printf("%.4f\n", m_b[i]);
//    }
    bool result = m_eq.solve(m_b, m_gamma);
//    cout << result << endl;
//    for (int i = 0; i < particles->size(); i++) {
//        printf("%.4f\n", m_gamma[i]);
//    }
//    cout << endl;
    temp.multiply(m_dp, m_gamma, particles->size() * 2, 1);

    for (int i = 0; i < particles->size(); i++) {
        Particle *p = particles->at(i);
        int n = m_counts[i];
        double mult = n > 0 ? (RELAXATION_PARAMETER / (double)n) : 0.,
               dx = m_dp[2 * i] * mult,
               dy = m_dp[2 * i + 1] * mult;
//        cout << dx << " " << dy << endl;

        p->ep.x += (fabs(dx) > EPSILON ? dx : 0);
        p->ep.y += (fabs(dy) > EPSILON ? dy : 0);

        if (stabile) {
            p->p.x += (fabs(dx) > EPSILON ? dx : 0);
            p->p.y += (fabs(dy) > EPSILON ? dy : 0);
        }
    }
}
