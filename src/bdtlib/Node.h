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

        float getErrorReduction();
        void setErrorReduction(float sErrorReduction);

        Node * getLeftDaughter();
        void setLeftDaughter(Node *sLeftDaughter);

        Node * getRightDaughter();
        void setRightDaughter(Node *sLeftDaughter);

        Node * getParent();
        void setParent(Node *sParent);

        float getSplitValue();
        void setSplitValue(float sSplitValue);

        int getSplitVariable();
        void setSplitVariable(int sSplitVar);

        float getFitValue();
        void setFitValue(float sFitValue);

        float getTotalError();
        void setTotalError(float sTotalError);

        float getAvgError();
        void setAvgError(float sAvgError);

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

        float splitValue;
        int splitVariable;

        float errorReduction;
        float totalError;
        float avgError;

        float fitValue;
        int numEvents;

        std::vector< std::vector<Event*> > events;
};

#endif
