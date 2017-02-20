// Node.h

#ifndef ADD_NODE
#define ADD_NODE

#include <string>
#include <vector>
#include "Event.h"
#include "SignificanceMetrics.h"

class Node
{
    public:
        Node();
        Node(std::string cName);
        ~Node();

        std::string getName();
        void setName(std::string sName);

        double getSignificanceGain();
        void setSignificanceGain(double sSignificanceGain);

        Node * getLeftDaughter();
        void setLeftDaughter(Node *sLeftDaughter);

        Node * getRightDaughter();
        void setRightDaughter(Node *sLeftDaughter);

        Node * getParent();
        void setParent(Node *sParent);

        double getSplitValue();
        void setSplitValue(double sSplitValue);

        int getSplitVariable();
        void setSplitVariable(int sSplitVar);

        double getSignificanceSquared();
        void setSignificanceSquared(double sSignificanceSquared);

        int getNumEvents();
        void setNumEvents(int sNumEvents);

        double getTotalSignal();
        void setTotalSignal(int sTotalSignal);

        double getTotalBackground();
        void setTotalBackground(int sTotalBackground);

        long long int getNumSignal();
        void setNumSignal(int sTotalSignal);

        long long int getNumBackground();
        void setNumBackground(int sTotalBackground);

        std::vector< std::vector<Event*> >& getEvents();
        void setEvents(std::vector< std::vector<Event*> >& sEvents);

        void calcOptimumSplit(SignificanceMetric* smetric, int nbins);
        void filterEventsToDaughters();
        Node* filterEventToDaughter(Event* e);
        void listEvents();
        void theMiracleOfChildBirth();
 
    private:
	std::string name;

        Node *leftDaughter;
        Node *rightDaughter;
        Node *parent;

        double splitValue;
        int splitVariable;

        double significanceGain;
        double significanceSquared;

        // sum of weights
        std::vector<double> totalSignalVec;     
        double totalSignal;

        std::vector<double> totalBackgroundVec;
        double totalBackground;

        // sum of signal/background training events used
        std::vector<long long int> numSignalVec;
        long long int numSignal;

        std::vector<long long int> numBackgroundVec;
        long long int numBackground;

        int numEvents;

        std::vector< std::vector<Event*> > events;
};

#endif
