#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <sstream>
#include <filesystem>

class SimpleProcess {
public:
    SimpleProcess(int isMC, int era, const std::string& sample) 
        : gensumweight(1.0), era(era), isMC(isMC), sample(sample), sumw(0.0) {
        
       
        if (gDirectory->FindObject("ht_reco")) {
            gDirectory->Delete("ht_reco");
        }
        if (gDirectory->FindObject("ht_reco_triggered")) {
            gDirectory->Delete("ht_reco_triggered");
        }

        ht_reco = new TH1F("ht_reco", "H_T [GeV];H_T [GeV];Events", 500, 0, 1000);
        ht_reco_triggered = new TH1F("ht_reco_triggered", "Triggered H_T [GeV];H_T [GeV];Events", 500, 0, 1000);
    }

    void process(TTree* tree) {
       
        std::vector<float>* Jet_pt = nullptr;
        std::vector<float>* Jet_eta = nullptr;
        std::vector<std::string>* hltResultName = nullptr;
        std::vector<bool>* hltResult = nullptr; 
        unsigned int scouting_trig = 0;

        tree->SetBranchAddress("Jet_pt", &Jet_pt);
        tree->SetBranchAddress("Jet_eta", &Jet_eta);

        if (tree->GetBranch("hltResultName")) {
            tree->SetBranchAddress("hltResultName", &hltResultName);
        } else {
            std::cout << "Error: Branch 'hltResultName' not found in the tree." << std::endl;
        }

        if (tree->GetBranch("hltResult")) {
            tree->SetBranchAddress("hltResult", &hltResult);
        } else {
            std::cout << "Error: Branch 'hltResult' not found in the tree." << std::endl;
        }
        tree->SetBranchAddress("scouting_trig", &scouting_trig);
        Long64_t nEntries = tree->GetEntries();
        std::cout << nEntries << std::endl;
        for (Long64_t i = 0; i < nEntries; ++i) {
            tree->GetEntry(i);
            std::vector<float> Ak4Jets;
            for (size_t j = 0; j < Jet_pt->size(); ++j) {
                if ((*Jet_pt)[j] > 30. && std::abs((*Jet_eta)[j]) < 4.7) {
                    Ak4Jets.push_back((*Jet_pt)[j]);
                }
            }
            std::vector<float> den;
            for (size_t k = 0; k < hltResultName->size(); ++k) {
                //std::cout << "Checking trigger: " << hltResultName->at(k) << " with result: " << hltResult->at(k) << std::endl;
                if (hltResult->at(k) == 1 && hltResultName->at(k).find("DST_DoubleMu3_noVtx_CaloScouting") != std::string::npos) {
                    den = Ak4Jets;

                }
            }

        float sum_den = 0.0;
        for (const auto& value : den) {
            sum_den += value;
        }
        ht_reco->Fill(sum_den);
        if (scouting_trig == 1) {
            ht_reco_triggered->Fill(sum_den);
        }
    }
    }

    void postprocess() {
        std::cout << "Total sum of weights for " << sample << ": " << sumw << std::endl;
    }

    TH1F* getHTReco() { return ht_reco; }
    TH1F* getHTRecoTriggered() { return ht_reco_triggered; }

private:
    float gensumweight;
    int era;
    int isMC;
    std::string sample;
    float sumw;
    TH1F* ht_reco;
    TH1F* ht_reco_triggered;
};

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <root_file>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    SimpleProcess processor(0, 2018, "test");

    std::cout << "Processing file: " << filename << std::endl;
    TFile input_file(filename.c_str());
    TTree* tree = (TTree*)input_file.Get("mmtree/tree");
    processor.process(tree);
    input_file.Close();

    processor.postprocess();

    std::string base_filename = filename.substr(filename.find_last_of('/') + 1);
    std::string output_filename = base_filename.substr(0, base_filename.find_last_of('.')) + "_out.root";
    std::cout << "Creating output file: " << output_filename << std::endl;
    TFile output_file(output_filename.c_str(), "RECREATE");
    if (output_file.IsZombie()) {
        std::cerr << "Error: Could not create output file " << output_filename << std::endl;
        return 1;
    }
    processor.getHTReco()->Write();
    processor.getHTRecoTriggered()->Write();
    output_file.Close();

    return 0;
}