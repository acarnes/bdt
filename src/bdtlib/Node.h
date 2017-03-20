// Node.h

#ifndef ADD_NODE
#define ADD_NODE

#include <string>
#include <vector>
#include "Event.h"
#include "SignificanceMetrics.hxx"

class Node
{
    public:
        Node();
        Node(std::string cName);
        ~Node();

        // bunch of get/set methods
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

        // signal weights
        double getTotalSignal();
        void setTotalSignal(double sTotalSignal);

        std::vector<double> getTotalSignalVec();
        void setTotalSignalVec(std::vector<double>& sTotalSignalVec);

        // bkg weights
        double getTotalBackground();
        void setTotalBackground(double sTotalBackground);

        std::vector<double> getTotalBackgroundVec();
        void setTotalBackgroundVec(std::vector<double>& sTotalBackgroundVec);

        // bkg out weights
        double getTotalBackgroundOut();
        void setTotalBackgroundOut(double sTotalBackgroundOut);

        // data out weights
        double getTotalDataOut();
        void setTotalDataOut(double sTotalDataOut);

        // num signal
        long long int getNumSignal();
        void setNumSignal(long long int sNumSignal);

        std::vector<long long int> getNumSignalVec();
        void setNumSignalVec(std::vector<long long int>& sNumSignalVec);

        // num bkg
        long long int getNumBackground();
        void setNumBackground(long long int sNumBackground);

        std::vector<long long int> getNumBackgroundVec();
        void setNumBackgroundVec(std::vector<long long int>& sNumBackgroundVec);

        // num bkg out
        long long int getNumBackgroundOut();
        void setNumBackgroundOut(long long int sNumBackground);

        // num data out
        long long int getNumDataOut();
        void setNumDataOut(long long int sNumData);

        std::vector< std::vector<Event*> >& getEvents();
        void setEvents(std::vector< std::vector<Event*> >& sEvents);

        // work functions
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
        double totalBackgroundOut;

        double totalDataOut;

        // sum of signal/background training events used
        std::vector<long long int> numSignalVec;
        long long int numSignal;

        std::vector<long long int> numBackgroundVec;
        long long int numBackground;
        long long int numBackgroundOut;

        long long int numDataOut;

        int numEvents;

        std::vector< std::vector<Event*> > events;
};

#endif
