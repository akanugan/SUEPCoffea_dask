import ROOT
import numpy as np
import os
import argparse  # Import argparse

class SimpleProcess:
    def __init__(self, isMC: int, era: int, sample: str):
        self.gensumweight = 1.0
        self.era = era
        self.isMC = isMC
        self.sample = sample

        if ROOT.gDirectory.FindObject("ht_reco"):
            ROOT.gDirectory.Delete("ht_reco")
        if ROOT.gDirectory.FindObject("ht_reco_triggered"):
            ROOT.gDirectory.Delete("ht_reco_triggered")

        self.ht_reco = ROOT.TH1F("ht_reco", "H_T [GeV];H_T [GeV];Events", 500, 0, 1000)
        self.ht_reco_triggered = ROOT.TH1F("ht_reco_triggered", "Triggered H_T [GeV];H_T [GeV];Events", 500, 0, 1000)
        self.sumw = 0.0

    def process(self, event):
        if hasattr(event, 'genWeight'):
            self.gensumweight = event.genWeight
            self.sumw += self.gensumweight
        else:
            self.sumw = 1.0

        # Create Ak4Jets
        Ak4Jets = []
        for i in range(len(event.Jet_pt)):
            if event.Jet_pt[i] > 30. and abs(event.Jet_eta[i]) < 4.7:
                Ak4Jets.append(event.Jet_pt[i])

        if event.scouting_trig == 1:
            Jets_triggered = Ak4Jets 
        else:
            Jets_triggered = []

        # Fill histograms
        self.ht_reco.Fill(sum(Ak4Jets))
        self.ht_reco_triggered.Fill(sum(Jets_triggered))

    def postprocess(self):
        print(f"Total sum of weights for {self.sample}: {self.sumw}")
        return {
            "ht_reco": self.ht_reco,
            "ht_reco_triggered": self.ht_reco_triggered,
            "sumw": self.sumw
        }

# Example usage
def main(data_dir):
#def main():
    processor = SimpleProcess(isMC=0, era=2018, sample='test')

    for filename in os.listdir(data_dir):
        if filename.endswith(".root"): 
            file_path = os.path.join(data_dir, filename)
            try:
                input_file = ROOT.TFile(file_path)  # Attempt to open the ROOT file
            except OSError as e:
                print(f"Warning: {e}. Skipping file: {file_path}")  # Log the error and skip the file
                continue  # Skip to the next file
            print(f"Processing file: {file_path}")

            #input_file = ROOT.TFile(file_path)
            tree = input_file.Get("mmtree/tree")
            for event_count, event in enumerate(tree):
                processor.process(event)
            input_file.Close()

    # with open('test.txt', 'r') as f:
    #     file_paths = f.read().splitlines()

    # for file_path in file_paths:
    #     if file_path.endswith(".root"): 
    #         print(f"Processing file: {file_path}")

    #         input_file = ROOT.TFile(file_path)
    #         tree = input_file.Get("mmtree/tree")
    #         for event_count, event in enumerate(tree):
    #             processor.process(event)
    #             input_file.Close()

    results = processor.postprocess()

    output_file = ROOT.TFile(args.out_file, "RECREATE")
    results["ht_reco"].Write()
    results["ht_reco_triggered"].Write()
    output_file.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process ROOT files.')
    parser.add_argument('-d', '--data_dir', type=str, help='Directory containing the ROOT files')
    parser.add_argument('-o', '--out_file', type=str)
    args = parser.parse_args()
    main(args.data_dir)
    #main()