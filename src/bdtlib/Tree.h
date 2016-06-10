// Tree.h

#ifndef ADD_TREE
#define ADD_TREE

#include <list>
#include "Node.h"
#include "TXMLEngine.h"
#include "tinyxml2.h"

class Tree
{
    public:
        Tree();
        Tree(std::vector< std::vector<Event*> >& cEvents);
        ~Tree();

        void setRootNode(Node *sRootNode);
        Node* getRootNode();

        void setTerminalNodes(std::list<Node*>& sTNodes);
        std::list<Node*>& getTerminalNodes();

        int getNumTerminalNodes();

        void buildTree(int nodeLimit);
        void calcError();
        void filterEvents(std::vector<Event*>& tEvents);
        void filterEventsRecursive(Node* node);
        Node* filterEvent(Event* e);
        Node* filterEventRecursive(Node* node, Event* e);

        void saveToXML(const char* filename);
        void saveToXMLRecursive(TXMLEngine* xml, Node* node, XMLNodePointer_t np);
        void addXMLAttributes(TXMLEngine* xml, Node* node, XMLNodePointer_t np);

        void tsaveToXML(const char* filename);
        void tsaveToXMLRecursive(tinyxml2::XMLDocument* xmlDoc, Node* node, tinyxml2::XMLElement* np);
        void taddXMLAttributes(Node* node, tinyxml2::XMLElement* np);

        void loadFromXML(const char* filename);
        void loadFromXMLRecursive(TXMLEngine* xml, XMLNodePointer_t node, Node* tnode);

        void rankVariables(std::vector<float>& v);
        void rankVariablesRecursive(Node* node, std::vector<float>& v);

        void getSplitValues(std::vector< std::vector<float> >& v);
        void getSplitValuesRecursive(Node* node, std::vector< std::vector<float> >& v);

    private:
        Node *rootNode;
        std::list<Node*> terminalNodes;
        int numTerminalNodes;
        float rmsError;
};

#endif
