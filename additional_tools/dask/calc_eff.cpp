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

    bool selectByL1(const std::vector<float>* Jet_pt, const std::vector<float>* Jet_eta, const std::vector<float>* Jet_phi, const std::vector<float>* Jet_m) {
        bool cutSingleJet = false;
        int countDoubleJet = 0;
        for (size_t i = 0; i < Jet_pt->size(); ++i) {
            if ((*Jet_pt)[i] > 250) cutSingleJet = true;
            if ((*Jet_pt)[i] > 50 && std::abs((*Jet_eta)[i]) < 2.5) countDoubleJet++;
        }

        if (cutSingleJet && countDoubleJet > 1) {
            for (size_t i = 0; i < Jet_pt->size(); ++i) {
                for (size_t j = i + 1; j < Jet_pt->size(); ++j) {
                    if (std::abs((*Jet_eta)[i] - (*Jet_eta)[j]) < 1.5) {
                        double p1 = std::sqrt(std::pow((*Jet_m)[i], 2) + std::pow(std::cosh((*Jet_eta)[i]) * (*Jet_pt)[i], 2));
                        double p2 = std::sqrt(std::pow((*Jet_m)[j], 2) + std::pow(std::cosh((*Jet_eta)[j]) * (*Jet_pt)[j], 2));
                        double kinematic = std::sqrt(std::pow(p1 + p2, 2) - std::pow((*Jet_pt)[i] * std::sin((*Jet_phi)[i]) + (*Jet_pt)[j] * std::sin((*Jet_phi)[j]), 2)
                            - std::pow((*Jet_pt)[i] * std::cos((*Jet_phi)[i]) + (*Jet_pt)[j] * std::cos((*Jet_phi)[j]), 2)
                            - std::pow(std::sinh((*Jet_eta)[i]) * (*Jet_pt)[i] + std::sinh((*Jet_eta)[j]) * (*Jet_pt)[j], 2));
                        if (kinematic > 420) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    void process(TTree* tree) {
       
        std::vector<Float_t>* Jet_pt = nullptr;
        std::vector<Float_t>* Jet_eta = nullptr;
        std::vector<float>* Jet_phi = nullptr;
        std::vector<float>* Jet_m = nullptr;
        std::vector<std::string>* hltResultName = nullptr;
        std::vector<bool>* hltResult = nullptr; 
        unsigned int scouting_trig = 0;

        tree->SetBranchAddress("Jet_pt", &Jet_pt);
        tree->SetBranchAddress("Jet_eta", &Jet_eta);
        tree->SetBranchAddress("Jet_phi", &Jet_phi);
        tree->SetBranchAddress("Jet_m", &Jet_m);

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
            std::vector<float> den; 
            if (selectByL1(Jet_pt, Jet_eta, Jet_phi, Jet_m)) {
                std::vector<float> Ak4Jets;
                for (size_t j = 0; j < Jet_pt->size(); ++j) {
                    if ((*Jet_pt)[j] > 30. && std::abs((*Jet_eta)[j]) < 4.7) {
                        Ak4Jets.push_back((*Jet_pt)[j]);
                    }
                }
                for (size_t k = 0; k < hltResultName->size(); ++k) {
                    //std::cout << "Checking trigger: " << hltResultName->at(k) << " with result: " << hltResult->at(k) << std::endl;
                    if (hltResult->at(k) == 1 && hltResultName->at(k).find("DST_DoubleMu3_noVtx_CaloScouting") != std::string::npos) {
                        den = Ak4Jets;
                    }
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

int main() {
    //std::string path = "/data/submit/akanugan/SUEP/data_PFComm_2017_flatNtuples/2017E/";
    //std::string path = "/data/submit/akanugan/SUEP/flatNtuples_MV/GluGluToSUEP_HT400_T0p35_mS200.000_mPhi1.400_T0.350_modehadronic_TuneCP5_13TeV/";
    // data 18 A and B got deleted from /work so copied other set of files to /data
    //std::string path = "data18D_PFComm/";
    std::string path = "/data/submit/akanugan/SUEP/data_PFComm_2018_flatNtuples/data18D/";

    SimpleProcess processor(0, 2018, "test");

    std::vector<std::string> files;
    for (const auto& entry : std::filesystem::directory_iterator(path)) {
        if (entry.path().extension() == ".root") {
            files.push_back(entry.path().string());
        }
    }

    for (const auto& filename : files) {
        std::cout << "Processing file: " << filename << std::endl;

        TFile input_file(filename.c_str());
        TTree* tree = (TTree*)input_file.Get("mmtree/tree");
        //TTree* tree = (TTree*)input_file.Get("Events");
        processor.process(tree);
        input_file.Close();
    }

    processor.postprocess();

    TFile output_file("data_18D_L1sel.root", "RECREATE");
    processor.getHTReco()->Write();
    processor.getHTRecoTriggered()->Write();
    output_file.Close();

    return 0;
}