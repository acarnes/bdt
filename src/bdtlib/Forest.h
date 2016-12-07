// Forest.h

#ifndef ADD_FOREST
#define ADD_FOREST

#include "Tree.h"
#include "LossFunctions.h"

class Forest
{
    public:

        // Constructor(s)/Destructor
        Forest();
        Forest(std::vector<Event*>& trainingEvents);
        ~Forest();

        // Get/Set
        void setTrainingEvents(std::vector<Event*>& trainingEvents);
        std::vector<Event*> getTrainingEvents();

        // Returns the number of trees in the forest.
        unsigned int size();

        // Get info on variable importance.
        void rankVariables(std::vector<int>& rank);
 
        // Output the list of split values used for each variable.
        void saveSplitValues(const char* savefilename);

        // Helpful operations
        void listEvents(std::vector< std::vector<Event*> >& e);
        void sortEventVectors(std::vector< std::vector<Event*> >& e);
        void loadForestFromXML(const char* directory, int numTrees); 

        // Perform the regression
        void fit(Tree *tree, double learningRate, LossFunction* l);
        void doRegression(int nodeLimit, int treeLimit, double learningRate, LossFunction* l, 
                          const char* savetreesdirectory, bool saveTrees);

        // Predict some events
        void updateEvents(Tree* tree);
        void appendCorrection(std::vector<Event*>& eventsp, int treenum);
        void predictEvents(std::vector<Event*>& eventsp, int trees);
        void appendCorrection(Event* e, int treenum);
        void predictEvent(Event* e, int trees);

        Tree* getTree(unsigned int i);

    private:

        std::vector< std::vector<Event*> > events;
        std::vector<Tree*> trees;
};

#endif
