// Node.h

#ifndef ADD_NODE
#define ADD_NODE

#include <string>
#include <vector>
#include "Event.h"

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

        std::vector< std::vector<Event*> >& getEvents();
        void setEvents(std::vector< std::vector<Event*> >& sEvents);

        void calcOptimumSplit();
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

        int numEvents;

        std::vector< std::vector<Event*> > events;
};

#endif
