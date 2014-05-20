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
        Forest(std::vector<Event*>& trainingEvents, std::vector<Event*>& testEvents);
        ~Forest();

        // Get/Set
        void setTrainingEvents(std::vector<Event*>& trainingEvents);
        void setTestEvents(std::vector<Event*>& testingEvents);
        std::vector<Event*> getTrainingEvents();
        std::vector<Event*> getTestEvents();

        // Info
        unsigned int size();


        // ------------------------------------------------------------------------------------------------
        // To Be Separated Out
        // ------------------------------------------------------------------------------------------------
        void readInTestingAndTrainingEvents(const char* inputfilename, LossFunction* l);
        void saveTestEvents(const char* savefilename);
        void loadUnknownEvents(const char* loadfile);
        void saveUnknownEvents(const char* savefilename);
        // ------------------------------------------------------------------------------------------------
        //
        // ------------------------------------------------------------------------------------------------


        // Helpful operations
        void listEvents(std::vector< std::vector<Event*> >& e);
        void sortEventVectors(std::vector< std::vector<Event*> >& e);
        void generate(Int_t numTrainEvents, Int_t numTestEvents, Double_t sigma);
        void loadForestFromXML(const char* directory, int numTrees); 

        // Perform the regression
        void doRegression(Int_t nodeLimit, Int_t treeLimit, Double_t learningRate, LossFunction* l, const char* savetreesdirectory, bool saveTrees, bool trackError, bool isTwoJets);
        void updateRegTargets(Tree *tree, Double_t learningRate, LossFunction* l);
        void doStochasticRegression(Int_t nodeLimit, Int_t treeLimit, Double_t learningRate, Double_t fraction, LossFunction* l);
        void prepareRandomSubsample(Double_t fraction);

        // Get info on variable importance.
        void rankVariables();

        // Predict testEvents/unknownEvents
        void resetTestEventPredictions();
        void predictTestEvents(Int_t trees);
        void predictUnknownEvents();
        void updateTestEvents(Tree* tree);
        void updateUnknownEvents(Tree* tree);

        // Calculate different types of error for a set of events(e.g. testEvents, unknownEvents, trainingEvents).
        Double_t returnRMS(std::vector<Event*>& v);
        Double_t returnResolution(std::vector<Event*>& v);
        Double_t returnDiscretizedAbsResolution(std::vector<Event*>& v, std::vector<Double_t> bins);

        // Keep track of different types of error as the tree is built.
        std::vector<Double_t> trainResolution;
        std::vector<Double_t> testResolution;
        std::vector<Double_t> trainRMS;
        std::vector<Double_t> testRMS;

    private:

        std::vector< std::vector<Event*> > events;
        std::vector< std::vector<Event*> > subSample;

        std::vector<Event*> testEvents;
        std::vector<Event*> unknownEvents;
        
        std::vector<Tree*> trees;
};

#endif
