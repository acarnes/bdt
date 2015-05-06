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
        void generate(Int_t numTrainEvents, Int_t numTestEvents, Double_t sigma);
        void loadForestFromXML(const char* directory, int numTrees); 

        // Perform the regression
        void updateRegTargets(Tree *tree, Double_t learningRate, LossFunction* l);
        void doRegression(Int_t nodeLimit, Int_t treeLimit, Double_t learningRate, LossFunction* l, 
                          const char* savetreesdirectory, bool saveTrees);

        // Stochastic Gradient Boosting
        void prepareRandomSubsample(Double_t fraction);
        void doStochasticRegression(Int_t nodeLimit, Int_t treeLimit, Double_t learningRate, 
                                    Double_t fraction, LossFunction* l);

        // Predict some events
        void updateEvents(Tree* tree);
        void appendCorrection(std::vector<Event*> eventsp, Int_t treenum);
        void predictEvents(std::vector<Event*> eventsp, Int_t trees);

    private:

        std::vector< std::vector<Event*> > events;
        std::vector< std::vector<Event*> > subSample;
        std::vector<Tree*> trees;
};

#endif
