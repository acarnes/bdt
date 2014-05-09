// Forest.h

#ifndef ADD_FOREST
#define ADD_FOREST

#include "Tree.h"
#include "LossFunctions.h"

class Forest
{
    public:
        Forest();
        ~Forest();

        void listEvents(std::vector< std::vector<Event*> >& e);
        void sortEventVectors(std::vector< std::vector<Event*> >& e);

        void generate(Int_t numTrainEvents, Int_t numTestEvents, Double_t sigma);
        void readInTestingAndTrainingEvents(const char* inputfilename, LossFunction* l);
        void readInTestingAndTrainingEventsFromDat(const char* inputfilename, LossFunction* l, bool isLog);
        void loadUnknownEvents(const char* loadfile);
        void loadForestFromXML(const char* directory, int numTrees); 

        void doRegression(Int_t nodeLimit, Int_t treeLimit, Double_t learningRate, LossFunction* l, const char* savetreesdirectory);

        void predictTestEvents();
        void predictUnknownEvents();

        void saveTestEvents(const char* savefilename);
        void saveTestEventsForJamie(const char* savefilename, bool isLog);
        void saveUnknownEvents(const char* savefilename);

        void doStochasticRegression(Int_t nodeLimit, Int_t treeLimit, Double_t learningRate, Double_t fraction, LossFunction* l);
        void prepareRandomSubsample(Double_t fraction);

        void updateTestEvents(Tree* tree);
        void updateUnknownEvents(Tree* tree);
        void updateRegTargets(Tree *tree, Double_t learningRate, LossFunction* l);

        Double_t returnRMS(std::vector<Event*>& v);
        Double_t returnResolution(std::vector<Event*>& v);
        Double_t returnInvResolution(std::vector<Event*>& v);
        void rankVariables();


    private:
        std::vector< std::vector<Event*> > events;
        std::vector< std::vector<Event*> > subSample;

        std::vector<Event*> testEvents;
        std::vector<Event*> unknownEvents;
        
        std::vector<Tree*> trees;
};

#endif
