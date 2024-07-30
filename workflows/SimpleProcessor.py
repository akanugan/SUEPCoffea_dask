

import awkward as ak
import numpy as np
import pandas as pd
import vector
vector.register_awkward()
from coffea import processor

class Simple_Process(processor.ProcessorABC):
    def __init__(self, isMC: int, era: int, sample: str) -> None:
        self.gensumweight = 1.0
        self.era = era
        self.isMC = isMC
        self.sample = sample

        self._accumulator = processor.dict_accumulator(
            {
                "ht_reco": hist.Hist(
                    "Events",
                    hist.Cat("dataset", "Dataset"),
                    hist.Bin("ht_reco", r"$H_T$ [GeV]", 50,0,2500),
                ),
                "ht_reco_triggered": hist.Hist(
                    "Events",
                    hist.Cat("dataset", "Dataset"),
                    hist.Bin("ht_reco_triggered", r"$H_T$ [GeV]",  50,0,2500),
                ),
                # "nmuons": hist.Hist(
                #     "Events",
                #     hist.Cat("dataset", "Dataset"),
                #     hist.Bin("nmuons", r"$N_{muons}$", 30, 0, 30),
                # ),
                # "muon_pt": hist.Hist(
                #     "Events",
                #     hist.Cat("dataset", "Dataset"),
                #     hist.Bin("muon_pt", r"$Muon p_{T}$ [GeV]", 10, 0, 200),
                # ),
                # "muon_pt_triggered": hist.Hist(
                #     "Events",
                #     hist.Cat("dataset", "Dataset"),
                #     hist.Bin("muon_pt_triggered", r"$Muon p_{T}$ [GeV]", 10, 0, 200),
                # ),
                # "MET": hist.Hist(
                #     "Events",
                #     hist.Cat("dataset", "Dataset"),
                #     hist.Bin("MET", r"$p_{T}^{miss}$ [GeV]", 50, 0, 200),
                # ),
                "sumw": processor.defaultdict_accumulator(float),
            }
        )
        
    @property
    def accumulator(self):
        return self._accumulator
    
    def process(self, events):
        output = self.accumulator
        dataset = events.metadata['dataset']
        #if "Muon" not in dataset:
        try:
            self.gensumweight = ak.sum(events.genWeight)
            output["sumw"][dataset] += ak.sum(events.genWeight)
        except:
            output["sumw"][dataset] = 1.0
        #output["sumw"][dataset] = 1.0
        
        # muons = ak.zip({
        #     "pt": events.Muon.pt,
        #     "eta": events.Muon.eta,
        #     "phi": events.Muon.phi,
        #     "mass": events.Muon.mass,
        #     "mediumId": events.Muon.mediumId
        # }, with_name="Momentum4D") 
        # muon_triggered = muons[events.HLT.IsoMu24 == 1]
        # muon_cut = (events.Muon.pt > 10) & \
        #     (abs(events.Muon.eta) <= 2.4) & \
        #     (events.Muon.mediumId == 1) 
        # muons = muons[muon_cut]
        # muon_triggered_cut = (muon_triggered.pt > 10) & \
        #     (abs(muon_triggered.eta) <= 2.4) & \
        #     (muon_triggered.mediumId == 1) 
        # muon_triggered = muon_triggered[muon_triggered_cut]
        
        Ak4Jets = ak.zip({
            "pt": events.Jet.pt,
            "eta": events.Jet.eta,
            "phi": events.Jet.phi,
            "mass": events.Jet.mass,
        }, with_name="Momentum4D")
        #Ak4Jets = Ak4Jets[events.HLT.Mu45_eta2p1 == 1]#for 2016
        
        if self.era == 2016:
            Jets_triggered = Ak4Jets[(events.HLT.PFHT900 == 1) & (events.HLT.Mu50 == 1)]
        else:
            #Jets_triggered = Ak4Jets[(events.HLT.PFHT1050 == 1) & (events.HLT.Mu50 == 1)]
            Jets_triggered = Ak4Jets[(events.scouting.trig == 1)] # scouting
        #Ak4Jets = Ak4Jets[events.HLT.Mu50 == 1]
        #Ak4Jets = Ak4Jets # all AK4 jets?
        Ak4JetCut = (Ak4Jets.pt > 30) & (abs(Ak4Jets.eta)<4.7)
        Ak4Jets = Ak4Jets[Ak4JetCut]
        Jets_triggeredCut = (Jets_triggered.pt > 30) & (abs(Jets_triggered.eta)<4.7)
        Jets_triggered = Jets_triggered[Jets_triggeredCut]
                
        # fill out hists

        output['ht_reco'].fill(ht_reco=ak.sum(Ak4Jets.pt,axis=-1), dataset=dataset)
        jet_trig = ak.to_numpy(ak.sum(Jets_triggered.pt,axis=-1),allow_missing=True)
        output['ht_reco_triggered'].fill(ht_reco_triggered=jet_trig, dataset=dataset)
        
        # output['nmuons'].fill(nmuons=ak.num(muons, axis=-1), dataset=dataset)
        # muons = muons[ak.num(muons, axis=-1)>0]
        # output['muon_pt'].fill(muon_pt=ak.max(muons.pt, axis=-1), dataset=dataset)
        # muon_trig = ak.to_numpy(ak.max(muon_triggered.pt, axis=-1),allow_missing=True)
        # output['muon_pt_triggered'].fill(muon_pt_triggered=muon_trig, dataset=dataset)
        # output['MET'].fill(MET=events.MET.pt, dataset=dataset)
                
        return output
        
    def postprocess(self, accumulator):
        return accumulator