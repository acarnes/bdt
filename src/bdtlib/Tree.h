// Tree.h

#ifndef ADD_TREE
#define ADD_TREE

#include <list>
#include "Node.h"
#include "tinyxml2.h"

class Tree
{
    public:
        Tree();
        Tree(std::vector<Event*>& cEvents);
        Tree(std::vector<Event*>& cEvents, int cnbins);
        ~Tree();

        void setTrainingEvents(std::vector<Event*>& trainingEvents);
        std::vector<Event*> getTrainingEvents();

        void setRootNode(Node *sRootNode);
        Node* getRootNode();

        void setTerminalNodes(std::list<Node*>& sTNodes);
        std::list<Node*>& getTerminalNodes();

        int getNumTerminalNodes();

        std::vector<std::string> getFeatureNames();
        void setFeatureNames(std::vector<std::string>& featureNames);

        void buildTree(int nodeLimit, SignificanceMetric* smetric);
        void calcSignificance();
        void filterEvents(std::vector<Event*>& tEvents);
        void filterEventsRecursive(Node* node);
        Node* filterEvent(Event* e);
        Node* filterEventRecursive(Node* node, Event* e);

        void saveToXML(const char* filename);
        void saveToXMLRecursive(tinyxml2::XMLDocument* xmlDoc, Node* node, tinyxml2::XMLElement* np);
        void addXMLAttributes(Node* node, tinyxml2::XMLElement* np);

        void loadFromXML(const char* filename);
        void loadFromXMLRecursive(tinyxml2::XMLElement* xmlnode, Node* tnode);

        void rankVariables(std::vector<double>& v);
        void rankVariablesRecursive(Node* node, std::vector<double>& v);
        void outputVariableRanking(std::vector<std::string>& rank);

        void getSplitValues(std::vector< std::vector<double> >& v);
        void getSplitValuesRecursive(Node* node, std::vector< std::vector<double> >& v);

        void sortEventVectors(std::vector< std::vector<Event*> >& e);
        void saveSplitValues(const char* savefilename);
        void listEvents(std::vector< std::vector<Event*> >& e);


    private:
        std::vector<std::string> featureNames;
        std::vector< std::vector<Event*> > events;
        Node *rootNode;
        std::list<Node*> terminalNodes;
        int numTerminalNodes;
        int nbins;
        double significance;
};

#endif
