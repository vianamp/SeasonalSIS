// ==============================================================
// SSIS: Seasonal SIS. Implementation of a time-dependent epidemic
// spreading model based on continuous time dynamics.
// ==============================================================

#include <list>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <algorithm>
#include <igraph/igraph.h>

double   rand_dbl() {
    return ((double)(rand()))/RAND_MAX;
}

long int rand_int(long int n) {
  long int rnd;
  do {
    rnd = rand();
  } while (rnd >= RAND_MAX);
  return rnd % n;
}

// Data structure to control how the transmissibility changes in time
//
// PRIVATE PARAMETERS:
// _t1:         1st piecewise time interval where the transimissibility increases
// _t2-_t1:     2nd piecewise time interval where the transimissibility decreases
// _lambda:     basal transimissibility
// _d_lambda:   total amount added to the transmissibility after _t1
//
// F1, F2:      piecewise functions that represent the transimissibility as function of time
// G1, G2:      integral of functions F1 and F2
// IG1, IG2     inverse of functions G1 and G2

class _transmissibility{
    private:
        double _t1, _t2;
        double _lambda;
        double _d_lambda;
        double _Lt1, _Lt2;

        double F1(double t);
        double F2(double t);

        double G1(double t);
        double G2(double t);

        double IG1(double t);
        double IG2(double t);

    public:
        ~_transmissibility() {}
        void Initialize(double t1, double t2, double lambda, double d_lambda);
        double Evaluate(double t);
        double EvaluateIntegral(double t);
        double EvaluateIntegralInverse(double L);
        void Print();

};

void _transmissibility::Initialize(double t1, double t2, double lambda, double d_lambda) {
    _t1 = t1;
    _t2 = t2;
    _lambda = lambda;
    _d_lambda = d_lambda;
    _Lt1 = G1(_t1);
    _Lt2 = _Lt1 + G2(_t2-_t1);
}

double _transmissibility::Evaluate(double t) {
    long int period = (long int)(t/_t2);
    double dt = t - period * _t2;
    if (dt < _t1)
        return F1(dt);
    else
        return F2(dt-_t1);
}

double _transmissibility::EvaluateIntegral(double t) {
    long int period = (long int)(t/_t2);
    double dt = t - period * _t2;
    if (dt < _t1)
        return period * _Lt2 + G1(dt);
    else
        return period * _Lt2 + _Lt1 + G2(dt-_t1);
}

double _transmissibility::EvaluateIntegralInverse(double L) {
    long int period = (long int)(L/_Lt2);
    double dL = L - period * _Lt2;
    if (dL < _Lt1)
        return period * _t2 + IG1(dL);
    else
        return period * _t2 + _t1 + IG2(dL-_Lt1);
}

void _transmissibility::Print() {
    double t = 0.0;
    FILE *f = fopen("temp.txt","w");
    fprintf(f,"t\tl\tL\n");
    while (t < 10.0) {
        fprintf(f,"%1.3f\t%1.3f\t%1.3f\n",t,Evaluate(t),EvaluateIntegral(t));
        t += 0.01;
    }
    fclose(f);
}

/* ================================================================
   PIECEWISE FUNCTION-DEPENDENT ROUTINES
=================================================================*/


double _transmissibility::F1(double t) {
    return _lambda;
}

double _transmissibility::F2(double t) {
    return _lambda + _d_lambda;
}

double _transmissibility::G1(double t) {
    return _lambda * t;
}

double _transmissibility::G2(double t) {
    return (_lambda + _d_lambda) * t;
}

double _transmissibility::IG1(double L) {
    return L / _lambda;
}

double _transmissibility::IG2(double L) {
    return L / (_lambda + _d_lambda);
}

/* ================================================================
   SIS
=================================================================*/

// Data structure to control the disease spreading
//
// PRIVATE PARAMETERS:
// _recovery:           Recovery rate, given by a fix value.
// _transimissibility:  Transmissibility rate, given by structure above.
//
// The method ImplementNextEvent uses a fast implementation of Gillespie algorithm
// to look for the next event, either infection or recovery. The event is then
// implemented and the time is updated.

class _SIS{
    private:
        double _mu;
        _transmissibility transmissibility;

    public:
        ~_SIS() {}
        _SIS(double t1, double t2, double lambda, double d_lambda, double tau);

        void Reset(igraph_t Graph);
        void InfectNode(igraph_t Graph, long int id);
        void RecoverNode(igraph_t Graph, long int id);

        void InfectRandomNode(igraph_t Graph);
        void InfectRandomNodes(double f, igraph_t Graph);

        long int ImplementNextEvent(igraph_t Graph);

        void RunSingleTrial(double f, igraph_t Graph, double Tmax, FILE *file, const char Text[]);

        double GetAsymptoticNumberOfInfectedNodes(double f, igraph_t Graph, int ntrials, double Tmax);

};

void _SIS::InfectRandomNode(igraph_t Graph) {
    long int node;
    do {
        node = rand_int(igraph_vcount(&Graph));
    } while (VAN(&Graph,"Infected",node)==0);
    InfectNode(Graph,node);
}

void _SIS::InfectRandomNodes(double f, igraph_t Graph) {
    long int i, j, k;
    std::vector<long int> Id;
    for (i = 0; i < igraph_vcount(&Graph); i++) {
        Id.push_back(i);
    }
    for (i = igraph_vcount(&Graph) - 1; i > 0; i--) {
        j = rand_int(i + 1);
        k = Id[j];
        Id[j] = Id[i];
        Id[i] = k;
    }
    k = (f > 0) ? (long int)(f*igraph_vcount(&Graph)) : 1;
    for (i = k; i--;) {
        InfectNode(Graph,Id[i]);
    }
    Id.clear();
}

void _SIS::RunSingleTrial(double f, igraph_t Graph, double Tmax, FILE *file, const char Text[]) {
    double t;
    Reset(Graph);
    long int c = 0, ninfected;
    InfectRandomNodes(f,Graph);
    do {
        t = GAN(&Graph,"t");
        ninfected = ImplementNextEvent(Graph);
        if (!(c%50)) {
            printf("Time = %1.3f\n",t);
            fprintf(file,"%s\t%1.3f\t%1.5f\n",Text,t,((double)ninfected)/igraph_vcount(&Graph));
            c = 0;
        }
        c++;
    } while ( (t < Tmax) & (ninfected > 0) );
}

double _SIS::GetAsymptoticNumberOfInfectedNodes(double f, igraph_t Graph, int ntrials, double Tmax) {
    double v = 0.0;
    long int ninfected;
    for (int trial = 0 ; trial < ntrials; trial++) {
        Reset(Graph);
        ninfected = 1;
        InfectRandomNodes(f,Graph);
        while ( (GAN(&Graph,"t") < Tmax) & (ninfected > 0) ) {
            ninfected = ImplementNextEvent(Graph);
        }
        v += ninfected;
    }
    return v / ntrials;
}

_SIS::_SIS(double t1, double t2, double lambda, double d_lambda, double mu) {
    _mu = mu;
    transmissibility.Initialize(t1,t2,lambda,d_lambda);
}

void _SIS::InfectNode(igraph_t Graph, long int id) {
    SETVAN(&Graph,"Infected",id,1);
}

void _SIS::RecoverNode(igraph_t Graph, long int id) {
    SETVAN(&Graph,"Infected",id,0);
}

void _SIS::Reset(igraph_t Graph) {
    SETGAN(&Graph,"t",0.0);
    SETGAN(&Graph,"L",0.0);
    for (long int node = igraph_vcount(&Graph); node--;)
        RecoverNode(Graph,node);
}

struct _event{
    public:
        char type;
        double p;
        long int node;
};

long int _SIS::ImplementNextEvent(igraph_t Graph) {

    double dt;
    long int i, j, node;
    igraph_vector_t Neigh;
    igraph_vector_init(&Neigh,0);
    std::vector<_event> Events;

    long int ninfected = 0;

    for (node = 0; node < igraph_vcount(&Graph); node++) {
        if (VAN(&Graph,"Infected",node) == 1) {

            ninfected++;

            _event event = {'r',1.0,node};
            Events.push_back(event);

            igraph_vector_clear(&Neigh);
            igraph_neighbors(&Graph, &Neigh, node, IGRAPH_ALL);
            for (i = 0; i < igraph_vector_size(&Neigh); i++) {
                j = (long int)VECTOR(Neigh)[i];
                if (VAN(&Graph,"Infected",j) == 0) {

                    _event event = {'i',2.0/200.0,j};
                    Events.push_back(event);

                }
            }
        }
    }

    double dtmin = 1E8;
    double totalp = 0.0;
    std::vector<_event>::iterator event;
    for(event=Events.begin();event!=Events.end();++event) {
        dt = -log( rand_dbl() ) / (*event).p;
        dtmin = (dt < dtmin) ? dtmin = dt : dtmin;
        totalp += (*event).p;
        (*event).p = totalp;
    }

    double p = 0;
    event = Events.begin();
    double r = totalp * rand_dbl();
    while (p < r) {
        p = (*event).p;
        ++event;
    }
    --event;

    double t = GAN(&Graph,"t");
    double L = GAN(&Graph,"L");

    SETGAN(&Graph,"t",t+dtmin);

    if ((*event).type=='i') {
        //SETGAN(&Graph,"L",L+r);
        InfectNode(Graph,(*event).node);
        ninfected++;
    } else {
        //SETGAN(&Graph,"L",transmissibility.EvaluateIntegral(t+smallest_recovery_interval));
        RecoverNode(Graph,(*event).node);
        ninfected--;
    }

    // Clear memory
    Events.clear();
    igraph_vector_destroy(&Neigh);

    return(ninfected);
}

/* ================================================================
   MAIN ROUTINE
=================================================================*/

igraph_t CreateLattice(int Lx, int Ly) {
    igraph_t Graph;
    igraph_vector_t dimvector;
    igraph_vector_init(&dimvector, 2);
    VECTOR(dimvector)[0] = Lx;
    VECTOR(dimvector)[1] = Ly;
    igraph_lattice(&Graph, &dimvector, 0, IGRAPH_UNDIRECTED, 0, 0);

    igraph_vector_t temp;
    igraph_vector_init(&temp, igraph_vcount(&Graph));
    igraph_vector_fill(&temp,0);
    SETVANV(&Graph,"Infected",&temp);
    igraph_vector_destroy(&temp);

    return Graph;
}

igraph_t CreateRandomGraph(int N, double p) {
    igraph_t Graph;
    igraph_erdos_renyi_game(&Graph, IGRAPH_ERDOS_RENYI_GNP, N, p, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    return Graph;
}

igraph_t CreateKRegularGraph(int N, int k) {
    igraph_t Graph;
    //igraph_k_regular_game(&Graph,N,k,false,false);
    return Graph;
}

igraph_t CreateCompleteGraph(int N) {
    igraph_t Graph;
    igraph_full(&Graph,N,false,false);
    return Graph;
}

int main(int argc, char *argv[]) {     

    srand(getpid());

    igraph_i_set_attribute_table(&igraph_cattribute_table);

    //igraph_t Graph = CreateLattice(10,10);
    //igraph_t Graph = CreateRandomGraph(400,0.05);
    igraph_t Graph = CreateCompleteGraph(200);

    FILE *fo = fopen("temp.txt","w");
    fprintf(fo,"model\ttime\ti\n");
    _SIS SISc(10,20,2.0,0.0,10.0); // T1, T2, L, DL, Tau
    //_SIS SISo(10,20,2.0,6.0,2.0);
    for (int run=1;run--;) {
        SISc.RunSingleTrial(1.0,Graph,100,fo,"Cont");
        //SISo.RunSingleTrial(1.0,Graph,50,fo,"Osci");
    }
    fclose(fo);

    return 0;
}

        
