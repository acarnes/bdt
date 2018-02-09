// Forest.h

#ifndef ADD_FOREST
#define ADD_FOREST

#include "Tree.h"
#include "SignificanceMetrics.hxx"

class Forest
{
    public:

        // Constructor(s)/Destructor
        Forest();
        Forest(std::vector<Event*>& trainingEvents);
        Forest(std::vector<Event*>& trainingEvents, double fEvents, int nFeatures);
        ~Forest();

        // Get/Set
        void setTrainingEvents(std::vector<Event*>& trainingEvents);
        std::vector<Event*> getTrainingEvents();
        std::vector<std::string> getFeatureNames();
        void setFeatureNames(std::vector<std::string>& cFeatureNames);

        // Returns the number of trees in the forest.
        unsigned int size();

        // Get info on variable importance.
        void appendFeatureRankings(Tree& tree);
        void outputFeatureRankings();
 
        // Output the list of split values used for each variable.
        void saveSplitValues(const char* savefilename);

        // Helpful operations
        void listEvents(std::vector< std::vector<Event*> >& e);
        void sortEventVectors(std::vector< std::vector<Event*> >& e);
        void loadForestFromXML(const char* directory, int numTrees); 

        // Perform the regression
        void doRegression(int nodeLimit, int treeLimit, int nbins, SignificanceMetric* s, 
                          const char* savetreesdirectory, bool saveTrees);

        // Predict some events
        void updateEvents(Tree* tree, int numtrees);
        void appendCorrection(std::vector<Event*>& eventsp, int treenum, int numtrees);
        void predictEvents(std::vector<Event*>& eventsp, int numtrees);
        void appendCorrection(Event* e, int treenum, int numtrees);
        void predictEvent(Event* e, int numtrees);

        Tree* getTree(unsigned int i);

    private:
        double fEvents; // fraction of events to use in each tree
        int nFeatures;  // number of features to use in each tree
        int numTrees;
        std::vector<double> featureRankings;
        std::vector<std::string> featureNames;
        std::vector< std::vector<Event*> > events;
        std::vector<Tree*> trees;
};

#endif
