/*          #################################################
           #              🏃️🏃‍♂️️🏃‍♀️️💨️                       #
          #  Di-dijets macro : 2023 ICPP                  #
         #  Refer to the CMS paper arXiv:2206.09997      #
        #               By Pilusa                       #
       #################################################


*/

#include <iostream>
#include <cmath>
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TCanvas.h" // Added include for TCanvas

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

void DIDIJETS()
{
  // Load the Delphes library
  gSystem->Load("libDelphes");

  // Create a chain of root trees
  TChain chain("Delphes");
  chain.Add("/home/electronicslab/Downloads/MG5_aMC_v3_4_2/4j_NLO/ATLAS/Events/run_03_CMS_CARD/tag_1_delphes_events.root");

  // Create an object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  // Book Histograms
     //..........Eta.............
  
  TH1F *histetaj1 = new TH1F("histEtaj1", "Etaj1;Etaj1 [rad];Event Count", 100, -3.15, 3.15); // Adjust the range as needed
  TH1F *histetaj2 = new TH1F("histEtaj2", "Etaj2;Etaj2 [rad];Event Count", 100, -3.15, 3.15); // Adjust the range as needed
  TH1F *histetaj3 = new TH1F("histEtaj3", "Etaj3;Etaj3 [rad];Event Count", 100, -3.15, 3.15); // Adjust the range as needed
  TH1F *histetaj4 = new TH1F("histEtaj4", "Etaj4;Etaj4 [rad];Event Count", 100, -3.15, 3.15); // Adjust the range as needed
   
          //.......Phi.........
  
  TH1F *histpj1 = new TH1F("histPhij1", "Phij1;Phij1 [rad];Event Count", 100, -3.15, 3.15); // Adjust the range as needed
  TH1F *histpj2 = new TH1F("histPhij2", "Phij2;Phij2 [rad];Event Count", 100, -3.15, 3.15); // Adjust the range as needed
  TH1F *histpj3 = new TH1F("histPhij3", "Phij3;Phij3 [rad];Event Count", 100, -3.15, 3.15); // Adjust the range as needed
  TH1F *histpj4 = new TH1F("histPhij4", "Phij4;Phij4 [rad];Event Count", 100, -3.15, 3.15); // Adjust the range as needed
  
  
    
       //..........changeR1 & changeR2..............
  
  TH1F *histR1 = new TH1F("histR1", "R1;R1 ;Event Count", 100, 0, 20); // Adjust the range as needed
  TH1F *histR2 = new TH1F("histR2", "R2;R2 ;Event Count", 100, 0, 20); // Adjust the range as needed
  
       //..........changeR..............
  
  TH1F *histR = new TH1F("histR", "R;R ;Event Count", 100, 0, 20); // Adjust the range as needed
        
         //....invariant mass ......
  
  TH1F *histM4j = new TH1F("histMjjjj", "M4j2800;M4j ;Event Count", 50, 0, 5000); // Adjust the range as needed
  TH1F *histMj1 = new TH1F("histMj1", "Mj1;Mj1 ;Event Count", 100, 0, 2000); // Adjust the range as needed
  TH1F *histMj2 = new TH1F("histMj2", "Mj2;Mj2 ;Event Count", 100, 0, 2000); // Adjust the range as needed
  
  
    
                   //......Sum of Eta1&Eta2, Eta3&Eta4, and Sum of all........
                   
    TH1F *histSETA12 = new TH1F("histSETA12", "#eta1 ;#eta_{1} [rad];Event Count", 100, -5, 5); // Adjust the range as needed              
    TH1F *histSETA34 = new TH1F("histSETA34", "#eta2 ;#eta_{2} [rad];Event Count", 100, -5, 5); // Adjust the range as needed
    TH1F *histSETA = new TH1F("histSETA", "#Delta#eta ;#Delta#eta [rad];Event Count", 100, -5, 5); // Adjust the range as needed
  
  
                   //......asymetry........
                   
    TH1F *histAsy = new TH1F("histAsy", "Asymetry;Assymetry ;Event Count", 100, -1, 1); // Adjust the range as needed
  
  
  
               //.........acceptance........
  TH1F *histAlpha = new TH1F("histMAcc", "M_{Z} = 3.5 TeV, M_{W} = 1 TeV ; #alpha ;Acceptance", 12, 0.0780743565300286, 0.34194470924690185); // Adjust the range as needed
  
  // Event counters
  Int_t totalEvents = 0;
  Int_t passedBothCuts = 0;
  Int_t passedPtJetCount = 0;
  Int_t passedEtaJetCount = 0;
  Int_t passedRCount = 0;
  Int_t passedSETACount = 0;
  Int_t passedAsyCount = 0;

  // Loop over all events
  for (Long64_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from the specified event
    treeReader->ReadEntry(entry);
    Event *event = (Event *)branchEvent->At(0);

    // Access jets from the Jet branch
    TClonesArray *jets = branchJet;

    // Increment the total event counter
    totalEvents++;

    // Check the number of jets in each event
    if (jets->GetEntries() >= 4)
    {
      // Count jets that pass the pT cut
  
      for (Int_t iJet = 0; iJet < jets->GetEntries(); ++iJet)
      {
        Jet *jet = (Jet *)jets->At(iJet);

        // Apply pT cut
        if (jet->PT > 80.0)
          passedPtJetCount++;
      }

      // Count jets that pass the eta cut
     
      for (Int_t iJet = 0; iJet < jets->GetEntries(); ++iJet)
      {
        Jet *jet = (Jet *)jets->At(iJet);

        // Apply eta cut
        if (std::abs(jet->Eta) < 2.5)
          passedEtaJetCount++;
      }

      // Check if the event passes both pT and eta cuts
      if (passedPtJetCount >= 4 && passedEtaJetCount >= 4)
      {
        passedBothCuts++;

        // Now you can fill histograms or perform other analyses on the selected events.
        // You can also increment the passedInvariantMassCut counter if your invariant mass cut is applied.

        // Declaring Lorentz Vectors for jet 1,2,3,4
        TLorentzVector vecJet1;
        TLorentzVector vecJet2;
        TLorentzVector vecJet3;
        TLorentzVector vecJet4;

        // Access the first four jets
        Jet *jet1 = (Jet *)jets->At(0);
        Jet *jet2 = (Jet *)jets->At(1);
        Jet *jet3 = (Jet *)jets->At(2);
        Jet *jet4 = (Jet *)jets->At(3);

        // Assign Lorentz Vectors to 4-momenta
        vecJet1 = jet1->P4();
        vecJet2 = jet2->P4();
        vecJet3 = jet3->P4();
        vecJet4 = jet4->P4();

        // Calculate eta for each jet
        Double_t etaj1 = vecJet1.Eta();
        Double_t etaj2 = vecJet2.Eta();
        Double_t etaj3 = vecJet3.Eta();
        Double_t etaj4 = vecJet4.Eta();

        histetaj1->Fill(etaj1);
        histetaj2->Fill(etaj2);
        histetaj3->Fill(etaj3);
        histetaj4->Fill(etaj4);

        // Calculate phi for each jet
        Double_t phij1 = vecJet1.Phi();
        Double_t phij2 = vecJet2.Phi();
        Double_t phij3 = vecJet3.Phi();
        Double_t phij4 = vecJet4.Phi();

        histpj1->Fill(phij1);
        histpj2->Fill(phij2);
        histpj3->Fill(phij3);
        histpj4->Fill(phij4);
        
        
//......THIS IS WHERE we START TO MAKE COMB1, COMB2, COMB3 event-by-event instead of doing one combination for all events..........
  //......We basically try to see which combination is the best for a certain event, instead of taking the best combination for the overall events........
     //......This means each event will have it's own comination.........

              


      



       
  // Loop over all possible combinations of changeRi
    Double_t minChangeR = 0;

            // Loop over all possible combinations of changeRi for the current event
            for (int i = 1; i <= 6; ++i)
            {
                Double_t changeRi1 = 0.0;
                Double_t changeRi2 = 0.0;

                // Calculate changeRi based on the combination i
                switch (i)
                {
                case 1:
                    changeRi1 = std::sqrt(std::pow(etaj1 - etaj2, 2) + std::pow(phij1 - phij2, 2));
                    changeRi2 = std::sqrt(std::pow(etaj3 - etaj4, 2) + std::pow(phij3 - phij4, 2));
                    break;
                case 2:
                    changeRi1 = std::sqrt(std::pow(etaj3 - etaj4, 2) + std::pow(phij3 - phij4, 2));
                    changeRi2 = std::sqrt(std::pow(etaj1 - etaj2, 2) + std::pow(phij1 - phij2, 2));
                    break;
                case 3:
                    changeRi1 = std::sqrt(std::pow(etaj1 - etaj3, 2) + std::pow(phij1 - phij3, 2));
                    changeRi2 = std::sqrt(std::pow(etaj2 - etaj4, 2) + std::pow(phij2 - phij4, 2));
                    break;
                case 4:
                    changeRi1 = std::sqrt(std::pow(etaj2 - etaj4, 2) + std::pow(phij2 - phij4, 2));
                    changeRi2 = std::sqrt(std::pow(etaj1 - etaj3, 2) + std::pow(phij1 - phij3, 2));
                    break;
                case 5:
                    changeRi1 = std::sqrt(std::pow(etaj1 - etaj4, 2) + std::pow(phij1 - phij4, 2));
                    changeRi2 = std::sqrt(std::pow(etaj2 - etaj3, 2) + std::pow(phij2 - phij3, 2));
                    break;
                case 6:
                    changeRi1 = std::sqrt(std::pow(etaj2 - etaj3, 2) + std::pow(phij2 - phij3, 2));
                    changeRi2 = std::sqrt(std::pow(etaj1 - etaj4, 2) + std::pow(phij1 - phij4, 2));
                    break;
                }

                // Calculate changeR based on the formula
                Double_t changeR = std::abs((changeRi1 - 0.8)) + std::abs((changeRi2 - 0.8));

                // Update minChangeR if the current combination has a smaller value
                if (changeR < std::numeric_limits<Double_t>::max())
                {
                    minChangeR = changeR;
                }
            }

// Now minChangeR contains the minimum changeR value



        // Apply the R cut
        if (minChangeR < 2)
        {
          passedRCount++;

          // Calculate Sum of ETA
         
           Double_t SETA12 = etaj1 + etaj4 ;
             histSETA12->Fill(SETA12);
           Double_t SETA34 = etaj2 + etaj3 ;
               histSETA34->Fill(SETA34);
           Double_t SETA = std::abs((etaj1 + etaj4) - (etaj2 + etaj3));
              histSETA->Fill(SETA); 
              
          // Apply SETA cut
          if (SETA < 1.1)
            passedSETACount++;

          // Calculate the Invariant Mass
          Double_t Mj1 = (vecJet1 + vecJet4).M();
          Double_t Mj2 = (vecJet2 + vecJet3).M();
          Double_t M4j = (vecJet1 + vecJet4 + vecJet2 + vecJet3).M();
          Double_t asy = std::abs((Mj1 - Mj2)) / (Mj1 + Mj2);
          histAsy->Fill(asy);
          if (Mj1 < 2000) histMj1->Fill(Mj1); 
          if (Mj2 < 2000) histMj2->Fill(Mj2); 
          if (M4j < 4000) histM4j->Fill(M4j); 

          // Apply the asymmetry cut
          if (asy < 0.1)
           {
            passedAsyCount++;
	  
          // Calculate Alpha
          Double_t Alpha = (Mj1 + Mj2)/(2*M4j);
          
          histAlpha->Fill(Alpha);
}
        }
      }
    }
  }

  // Show resulting histograms
  
  
  /*
  
  //.........Plot Eta...............
  TCanvas *canvas = new TCanvas("Eta", "Eta");
  canvas->Divide(2, 2);
  canvas->cd(1);
  histetaj1->Draw();
    canvas->cd(2);
  histetaj2->Draw();
    canvas->cd(3);
  histetaj3->Draw();
    canvas->cd(4);
  histetaj4->Draw();
  
  */
  
  /*
    //.........Plot Phi...............
  TCanvas *canvas = new TCanvas("Phi", "Phi");
  canvas->Divide(2, 2);
  canvas->cd(1);
  histpj1->Draw();
    canvas->cd(2);
  histpj2->Draw();
    canvas->cd(3);
  histpj3->Draw();
    canvas->cd(4);
  histpj4->Draw();
  
*/  
  
  
  
 /* //.........Plot ChangeR1 and ChangeR2
  
    TCanvas *canvas = new TCanvas("changeR1&R2", "changeR1&R2");
  canvas->Divide(2, 2);
  canvas->cd(1);
  histR1->Draw();
    canvas->cd(2);
  histR2->Draw();
*/
  
  /*
    //.........Plot ChangeR for different jet pairs
  
    TCanvas *canvas = new TCanvas("changeR", "changeR");
   canvas->Divide(1, 1);
   canvas->cd(1);
   histR->Draw();
  */
 
    //.........Plot Sum of Eta1&Eta2, Eta3&Eta4, and Sum of all.........
  
   /* TCanvas *canvas = new TCanvas("Sum of Eta", "Sum of Eta");
   canvas->Divide(2, 3);
   canvas->cd(1);
   histSETA12->Draw();
   canvas->cd(2);
   histSETA34->Draw();
   canvas->cd(3);
   histSETA->Draw();
  */
  
 //.......plot invariant mass
 /*
   TCanvas *canvas = new TCanvas("M", "M");
   canvas->Divide(2, 3);
   canvas->cd(1);
    histMj1->Draw();
   canvas->cd(2);
    histMj2->Draw();
   canvas->cd(3);
    histM4j->Draw();
    
    */
    
/*
      //.........Plot asymmetry
  
    TCanvas *canvas = new TCanvas("asymetry", "asymetry");
   canvas->Divide(1, 1);
   canvas->cd(1);
   histAsy->Draw(); 
   
 */  
   

  
//            >>>>>>>>>>>> This part is for the acceptance >>>>>>>>>

   TCanvas *canvas = new TCanvas("Acc", "angular distn");
  canvas->Divide(1, 1);
  canvas->cd(1);
  histAlpha->GetYaxis()->SetRangeUser(0, histAlpha->GetMaximum() + histAlpha->GetMaximum() * 0.4); // Adjust the range as needed
  histAlpha->GetYaxis()->SetTitle("Acceptance"); // Set the y-axis title
  histAlpha->GetYaxis()->CenterTitle(); // Center the y-axis title
  histAlpha->Scale(1.0 / 100000); // Scale the histogram values
  histAlpha->Draw();

  // Print the number of events at each stage
  std::cout << "Total events: " << totalEvents << std::endl;
  std::cout << "Events passing pT cut: " << passedPtJetCount << std::endl;
  std::cout << "Events passing Eta cut: " << passedEtaJetCount << std::endl;
  std::cout << "Events passing both pT and eta cuts: " << passedBothCuts << std::endl;
  std::cout << "Events passing R < 2: " << passedRCount << std::endl;
  std::cout << "Events passing Sum of ETA < 1.1: " << passedSETACount << std::endl;
  std::cout << "Events passing Asymetry < 0.1: " << passedAsyCount << std::endl;
}

