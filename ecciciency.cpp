#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"

std::vector<double> CreateLogBinning(int nbins, double xmin, double xmax) {
    std::vector<double> bin_edges(nbins + 1);
    double logxmin = std::log10(xmin);
    double logxmax = std::log10(xmax);
    double bin_width = (logxmax - logxmin) / nbins;
    for (int i = 0; i <= nbins; ++i) {
        bin_edges[i] = std::pow(10, logxmin + i * bin_width);
    }
    return bin_edges;
}

int gatto() {
    const char* combinedFileName = "A.combinedhistAbs.root";

    // Apri il file combinato per analisi
    TFile *combinedFile = TFile::Open(combinedFileName, "UPDATE");
    if (!combinedFile || combinedFile->IsZombie()) {
        std::cerr << "Error: Cannot open combined file " << combinedFileName << std::endl;
        return 1;
    }

    /*
    // Creare la cartella per salvare le immagini
    const char* imgDir = "canvas_Eff";
    if (mkdir(imgDir, 0777) && errno != EEXIST) {
        std::cerr << "Error: Cannot create directory " << imgDir << std::endl;
        combinedFile->Close();
        return 1;
    }
    */

    // _________________________________________________________________ Q2

    TH1D *pion11 = (TH1D*)combinedFile->Get("truepionQ22");
    TH1D *pion12 = (TH1D*)combinedFile->Get("RecpionQ22");
    TH1D *pion13 = (TH1D*)combinedFile->Get("RealpionQ2");
    TH1D *kaon11 = (TH1D*)combinedFile->Get("truekaonQ22");
    TH1D *kaon12 = (TH1D*)combinedFile->Get("ReckaonQ22");
    TH1D *kaon13 = (TH1D*)combinedFile->Get("RealkaonQ2");
    TH1D *proton11 = (TH1D*)combinedFile->Get("trueprotonQ22");
    TH1D *proton12 = (TH1D*)combinedFile->Get("RecprotonQ22");

    // __________________________________________________________________ x_bj

    TH1D *pion21 = (TH1D*)combinedFile->Get("truepion_xbj2");
    TH1D *pion22 = (TH1D*)combinedFile->Get("Recpion_xbj2");
    TH1D *pion23 = (TH1D*)combinedFile->Get("RealpionX");
    TH1D *kaon21 = (TH1D*)combinedFile->Get("truekaon_xbj2");
    TH1D *kaon22 = (TH1D*)combinedFile->Get("Reckaon_xbj2");
    TH1D *kaon23 = (TH1D*)combinedFile->Get("RealkaonX");
    TH1D *proton21 = (TH1D*)combinedFile->Get("trueproton_xbj2");
    TH1D *proton22 = (TH1D*)combinedFile->Get("Recproton_xbj2");

    // __________________________________________________________________ z

    TH1D *pion31 = (TH1D*)combinedFile->Get("truepion_z2");
    TH1D *pion32 = (TH1D*)combinedFile->Get("Recpion_z2");
    TH1D *pion33 = (TH1D*)combinedFile->Get("RealpionZ");
    TH1D *kaon31 = (TH1D*)combinedFile->Get("truekaon_z2");
    TH1D *kaon32 = (TH1D*)combinedFile->Get("Reckaon_z2");
    TH1D *kaon33 = (TH1D*)combinedFile->Get("RealkaonZ");
    TH1D *proton31 = (TH1D*)combinedFile->Get("trueproton_z2");
    TH1D *proton32 = (TH1D*)combinedFile->Get("Recproton_z2");

    // __________________________________________________________________ P_hT

    TH1D *pion41 = (TH1D*)combinedFile->Get("truepion_PhT2");
    TH1D *pion42 = (TH1D*)combinedFile->Get("Recpion_PhT2");
    TH1D *pion43 = (TH1D*)combinedFile->Get("RealpionPht");
    TH1D *kaon41 = (TH1D*)combinedFile->Get("truekaon_PhT2");
    TH1D *kaon42 = (TH1D*)combinedFile->Get("Reckaon_PhT2");
    TH1D *kaon43 = (TH1D*)combinedFile->Get("RealkaonPhT");
    TH1D *proton41 = (TH1D*)combinedFile->Get("trueproton_PhT2");
    TH1D *proton42 = (TH1D*)combinedFile->Get("Recproton_PhT2");

    // __________________________________________________________________ Eta

    TH1D *pion51 = (TH1D*)combinedFile->Get("truepionEta2");
    TH1D *pion52 = (TH1D*)combinedFile->Get("RecpionEta2");
    TH1D *pion53 = (TH1D*)combinedFile->Get("RealpionEta");
    TH1D *kaon51 = (TH1D*)combinedFile->Get("truekaonEta2");
    TH1D *kaon52 = (TH1D*)combinedFile->Get("ReckaonEta2");
    TH1D *kaon53 = (TH1D*)combinedFile->Get("RealkaonEta");
    TH1D *proton51 = (TH1D*)combinedFile->Get("trueprotonEta2");
    TH1D *proton52 = (TH1D*)combinedFile->Get("RecprotonEta2");

    // __________________________________________________________________ Phi

    TH1D *pion61 = (TH1D*)combinedFile->Get("truepionPhi2");
    TH1D *pion62 = (TH1D*)combinedFile->Get("RecpionPhi2");
    TH1D *pion63 = (TH1D*)combinedFile->Get("RealpionPhi");
    TH1D *kaon61 = (TH1D*)combinedFile->Get("truekaonPhi2");
    TH1D *kaon62 = (TH1D*)combinedFile->Get("ReckaonPhi2");
    TH1D *kaon63 = (TH1D*)combinedFile->Get("RealkaonPhi");
    TH1D *proton61 = (TH1D*)combinedFile->Get("trueprotonPhi2");
    TH1D *proton62 = (TH1D*)combinedFile->Get("RecprotonPhi2");

    // __________________________________________________________________ Mom

    TH1D *pion71 = (TH1D*)combinedFile->Get("truepion_mom2");
    TH1D *pion72 = (TH1D*)combinedFile->Get("Recpion_mom2");
    TH1D *pion73 = (TH1D*)combinedFile->Get("RealpionMom");
    TH1D *kaon71 = (TH1D*)combinedFile->Get("truekaon_mom2");
    TH1D *kaon72 = (TH1D*)combinedFile->Get("Reckaon_mom2");
    TH1D *kaon73 = (TH1D*)combinedFile->Get("RealkaonMom");
    TH1D *proton71 = (TH1D*)combinedFile->Get("trueproton_mom2");
    TH1D *proton72 = (TH1D*)combinedFile->Get("Recproton_mom2");

    // ________________________________________________________________________________________________________________________________

    // _________________________________________________________________ Q2

    TH1D *fpion11 = (TH1D*)combinedFile->Get("truepionQ22");
    TH1D *fpion12 = (TH1D*)combinedFile->Get("RecpionQ22");
    TH1D *fpion13 = (TH1D*)combinedFile->Get("RealpionQ2");
    TH1D *fkaon11 = (TH1D*)combinedFile->Get("truekaonQ22");
    TH1D *fkaon12 = (TH1D*)combinedFile->Get("ReckaonQ22");
    TH1D *fkaon13 = (TH1D*)combinedFile->Get("RealkaonQ2");
    TH1D *fproton11 = (TH1D*)combinedFile->Get("trueprotonQ22");
    TH1D *fproton12 = (TH1D*)combinedFile->Get("RecprotonQ22");

    // __________________________________________________________________ x_bj

    TH1D *fpion21 = (TH1D*)combinedFile->Get("truepion_xbj2");
    TH1D *fpion22 = (TH1D*)combinedFile->Get("Recpion_xbj2");
    TH1D *fpion23 = (TH1D*)combinedFile->Get("RealpionX");
    TH1D *fkaon21 = (TH1D*)combinedFile->Get("truekaon_xbj2");
    TH1D *fkaon22 = (TH1D*)combinedFile->Get("Reckaon_xbj2");
    TH1D *fkaon23 = (TH1D*)combinedFile->Get("RealkaonX");
    TH1D *fproton21 = (TH1D*)combinedFile->Get("trueproton_xbj2");
    TH1D *fproton22 = (TH1D*)combinedFile->Get("Recproton_xbj2");

    // __________________________________________________________________ z

    TH1D *fpion31 = (TH1D*)combinedFile->Get("truepion_z2");
    TH1D *fpion32 = (TH1D*)combinedFile->Get("Recpion_z2");
    TH1D *fpion33 = (TH1D*)combinedFile->Get("RealpionZ");
    TH1D *fkaon31 = (TH1D*)combinedFile->Get("truekaon_z2");
    TH1D *fkaon32 = (TH1D*)combinedFile->Get("Reckaon_z2");
    TH1D *fkaon33 = (TH1D*)combinedFile->Get("RealkaonZ");
    TH1D *fproton31 = (TH1D*)combinedFile->Get("trueproton_z2");
    TH1D *fproton32 = (TH1D*)combinedFile->Get("Recproton_z2");

    // __________________________________________________________________ P_hT

    TH1D *fpion41 = (TH1D*)combinedFile->Get("truepion_PhT2");
    TH1D *fpion42 = (TH1D*)combinedFile->Get("Recpion_PhT2");
    TH1D *fpion43 = (TH1D*)combinedFile->Get("RealpionPht");
    TH1D *fkaon41 = (TH1D*)combinedFile->Get("truekaon_PhT2");
    TH1D *fkaon42 = (TH1D*)combinedFile->Get("Reckaon_PhT2");
    TH1D *fkaon43 = (TH1D*)combinedFile->Get("RealkaonPhT");
    TH1D *fproton41 = (TH1D*)combinedFile->Get("trueproton_PhT2");
    TH1D *fproton42 = (TH1D*)combinedFile->Get("Recproton_PhT2");

    // __________________________________________________________________ Eta

    TH1D *fpion51 = (TH1D*)combinedFile->Get("truepionEta2");
    TH1D *fpion52 = (TH1D*)combinedFile->Get("RecpionEta2");
    TH1D *fpion53 = (TH1D*)combinedFile->Get("RealpionEta");
    TH1D *fkaon51 = (TH1D*)combinedFile->Get("truekaonEta2");
    TH1D *fkaon52 = (TH1D*)combinedFile->Get("ReckaonEta2");
    TH1D *fkaon53 = (TH1D*)combinedFile->Get("RealkaonEta");
    TH1D *fproton51 = (TH1D*)combinedFile->Get("trueprotonEta2");
    TH1D *fproton52 = (TH1D*)combinedFile->Get("RecprotonEta2");

    // __________________________________________________________________ Phi

    TH1D *fpion61 = (TH1D*)combinedFile->Get("truepionPhi2");
    TH1D *fpion62 = (TH1D*)combinedFile->Get("RecpionPhi2");
    TH1D *fpion63 = (TH1D*)combinedFile->Get("RealpionPhi");
    TH1D *fkaon61 = (TH1D*)combinedFile->Get("truekaonPhi2");
    TH1D *fkaon62 = (TH1D*)combinedFile->Get("ReckaonPhi2");
    TH1D *fkaon63 = (TH1D*)combinedFile->Get("RealkaonPhi");
    TH1D *fproton61 = (TH1D*)combinedFile->Get("trueprotonPhi2");
    TH1D *fproton62 = (TH1D*)combinedFile->Get("RecprotonPhi2");

    // __________________________________________________________________ Mom

    TH1D *fpion71 = (TH1D*)combinedFile->Get("truepion_mom2");
    TH1D *fpion72 = (TH1D*)combinedFile->Get("Recpion_mom2");
    TH1D *fpion73 = (TH1D*)combinedFile->Get("RealpionMom");
    TH1D *fkaon71 = (TH1D*)combinedFile->Get("truekaon_mom2");
    TH1D *fkaon72 = (TH1D*)combinedFile->Get("Reckaon_mom2");
    TH1D *fkaon73 = (TH1D*)combinedFile->Get("RealkaonMom");
    TH1D *fproton71 = (TH1D*)combinedFile->Get("trueproton_mom2");
    TH1D *fproton72 = (TH1D*)combinedFile->Get("Recproton_mom2");


    TFile *outputFile = TFile::Open("B.efficienciesAbs.root", "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Cannot create output file." << std::endl;
        combinedFile->Close();
        return 1;
    }


    // ______________________________________________________________ Q2 _____________________________________________________________

    int nbins = 80;
    double xmin_xbj = 5e-3;
    double xmax_xbj = 1;
    double xmin_Q2 = 0.9;
    double xmax_Q2 = 100.;
    std::vector<double> log_bins_Q2 = CreateLogBinning(nbins, xmin_Q2, xmax_Q2);
    std::vector<double> log_bins_xbj = CreateLogBinning(nbins, xmin_xbj, xmax_xbj);
    TH1D *pionQ2 = new TH1D("pionQ2", "Efficiency reconstruction with MC ID  |  Q^2  |  18x275 GeV", nbins, log_bins_Q2.data());
    pionQ2->Divide(pion12, pion11);
    pionQ2->SetLineColor(kRed);
    pionQ2->SetStats(kFALSE);
    pionQ2->GetXaxis()->SetTitle("Q^2 [GeV^2]");
    pionQ2->GetYaxis()->SetTitle("Rec / MC events");
    pionQ2->GetYaxis()->SetRangeUser(0, 1.2);
    TH1D *kaonQ2 = new TH1D("kaonQ2", "Ratio of Histograms", nbins, log_bins_Q2.data());
    kaonQ2->Divide(kaon12,kaon11);
    kaonQ2->SetLineColor(kBlue);
    TH1D *protonQ2 = new TH1D("protonQ2", "Ratio of Histograms", nbins, log_bins_Q2.data());
    protonQ2->Divide(proton12,proton11);
    protonQ2->SetLineColor(kOrange);

    TLegend *legend2 = new TLegend(0.75, 0.78, 0.88, 0.88);
    legend2->AddEntry(pionQ2, "Pions", "l");
    legend2->AddEntry(kaonQ2, "Kaons", "l");
    legend2->AddEntry(protonQ2, "Protons", "l");

    TCanvas *Ratio_Q2 = new TCanvas("Ratio_Q2", "Ratio Canvas", 800, 600);
    Ratio_Q2->SetLogx();
    pionQ2->Draw("HIST");
    kaonQ2->Draw("HIST SAME");
    protonQ2->Draw("HIST SAME");
    legend2->Draw();
    
    Ratio_Q2->Write();
    //Ratio_Q2->SaveAs("canvas_Eff/Ratio_Q2.jpg");

    TH1D *pionQ2_ID = new TH1D("pionQ2_ID", "Reconstruction x PID efficiency  |  Q^2  |  18x275 GeV", nbins, log_bins_Q2.data());
    pionQ2_ID->Divide(pion13, pion11);
    pionQ2_ID->SetLineColor(kRed);
    pionQ2_ID->SetLineWidth(1);
    pionQ2_ID->SetStats(kFALSE);
    pionQ2_ID->GetXaxis()->SetTitle("Q^2 [GeV^2]");
    pionQ2_ID->GetYaxis()->SetTitle("Rec / MC events");
    pionQ2_ID->GetYaxis()->SetRangeUser(0, 1.2);
    TH1D *kaonQ2_ID = new TH1D("kaonQ2_ID", "Reconstruction x PID efficiency  |  Q^2  |  18x275 GeV", nbins, log_bins_Q2.data());
    kaonQ2_ID->Divide(kaon13,kaon11);
    kaonQ2_ID->SetLineColor(kBlue);
    kaonQ2_ID->GetXaxis()->SetTitle("Q^2 [GeV^2]");
    kaonQ2_ID->GetYaxis()->SetTitle("Rec / MC events");
    kaonQ2_ID->GetYaxis()->SetRangeUser(0, 1.2);
    kaonQ2_ID->SetStats(kFALSE);

    TLegend *legend = new TLegend(0.75, 0.78, 0.88, 0.88);
    legend->AddEntry(pionQ2_ID, "Pions", "l");
    legend->AddEntry(kaonQ2_ID, "Kaons", "l");
    //legend2->AddEntry(protonQ2_ID, "Protons", "l");

    TCanvas *Ratio_Q2_ID = new TCanvas("Ratio_Q2_ID", "Ratio Canvas", 800, 600);
    Ratio_Q2_ID->SetLogx();
    kaonQ2_ID->Draw("HIST");
    pionQ2_ID->Draw("HIST SAME");
    //protonQ2_ID->Draw("HIST SAME");
    legend->Draw();
    
    Ratio_Q2_ID->Write();
    //Ratio_Q2_ID->SaveAs("canvas_Eff/Ratio_Q2.jpg");

    // FRACTION____________________________________________________________________________________________________
    TCanvas *MCfraction_Q2 = new TCanvas("MCfraction_Q2", "Production with Q2", 800, 600);
    MCfraction_Q2->SetLogx();

    int nBins = fpion11->GetNbinsX();
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = fpion11->GetBinContent(i) + fkaon11->GetBinContent(i) + fproton11->GetBinContent(i);

        if (sumBinContents > 0) {
            fpion11->SetBinContent(i, fpion11->GetBinContent(i) / (sumBinContents));
            fkaon11->SetBinContent(i, fkaon11->GetBinContent(i) / (sumBinContents));
            fproton11->SetBinContent(i, fproton11->GetBinContent(i) / (sumBinContents));
        }
    }
    fpion11->SetStats(kFALSE);
    fpion11->SetLineColor(kRed);
    fpion11->SetTitle("MC Production of charged particles | Q^2 | 18x275 GeV");
    fpion11->GetXaxis()->SetTitle("Q^2 [GeV^2]");
    fpion11->GetYaxis()->SetRangeUser(0, 1.);
    fpion11->GetYaxis()->SetTitle("Relative fraction");
    fpion11->Draw("HIST");
    fkaon11->SetLineColor(kBlue);
    fkaon11->Draw("HIST SAME");
    fproton11->SetLineColor(kOrange);
    fproton11->Draw("HIST SAME");
    legend2->Draw();

    //MCfraction_Q2->SetTitle("Production of charged particles with Q2");
    MCfraction_Q2->Update();
    MCfraction_Q2->Write();

    // RECONSTRUCTED FRACTION
    TCanvas *RecFraction_Q2 = new TCanvas("RecFraction_Q2", "Production with Q2", 800, 600);
    RecFraction_Q2->SetLogx();

    //int nBins = fpion12->GetNbinsX();
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = fpion12->GetBinContent(i) + fkaon12->GetBinContent(i) + fproton12->GetBinContent(i);

        if (sumBinContents > 0) {
            fpion12->SetBinContent(i, fpion12->GetBinContent(i) / (sumBinContents));
            fkaon12->SetBinContent(i, fkaon12->GetBinContent(i) / (sumBinContents));
            fproton12->SetBinContent(i, fproton12->GetBinContent(i) / (sumBinContents));
        }
    }
    fpion12->SetStats(kFALSE);
    fpion12->SetLineColor(kRed);
    fpion12->SetTitle("Reconstructed production of charged particles | Q^2 | 18x275 GeV");
    fpion12->GetXaxis()->SetTitle("Q^2 [GeV^2]");
    fpion12->GetYaxis()->SetRangeUser(0, 1.);
    fpion12->GetYaxis()->SetTitle("Relative fraction");
    fpion12->Draw("HIST");
    fkaon12->SetLineColor(kBlue);
    fkaon12->Draw("HIST SAME");
    fproton12->SetLineColor(kOrange);
    fproton12->Draw("HIST SAME");
    legend2->Draw();

    //RecFraction_Q2->SetTitle("Production of charged particles with Q2");
    RecFraction_Q2->Update();
    RecFraction_Q2->Write();


    // ______________________________________________________________ x_bj _____________________________________________________________

    TH1D *pion_xbj = new TH1D("pion_xbj", "Efficiency reconstruction with MC ID  |  x_B  |  18x275 GeV", nbins, log_bins_xbj.data());
    pion_xbj->Divide(pion22, pion21);
    pion_xbj->SetLineColor(kRed);
    pion_xbj->SetStats(kFALSE);
    pion_xbj->GetXaxis()->SetTitle("x_B");
    pion_xbj->GetYaxis()->SetTitle("Rec / MC events");
    pion_xbj->GetYaxis()->SetRangeUser(0, 1.2);
    TH1D *kaon_xbj = new TH1D("kaon_xbj", "Ratio of Histograms", nbins, log_bins_xbj.data());
    kaon_xbj->Divide(kaon22,kaon21);
    kaon_xbj->SetLineColor(kBlue);
    TH1D *proton_xbj = new TH1D("proton_xbj", "Ratio of Histograms", nbins, log_bins_xbj.data());
    proton_xbj->Divide(proton22,proton21);
    proton_xbj->SetLineColor(kOrange);

    TCanvas *Ratio_xbj = new TCanvas("Ratio_xbj", "Ratio Canvas", 800, 600);
    Ratio_xbj->SetLogx();
    pion_xbj->Draw("HIST");
    kaon_xbj->Draw("HIST SAME");
    proton_xbj->Draw("HIST SAME");
    legend2->Draw();
    
    Ratio_xbj->Write();

    TH1D *pion_xbj_ID = new TH1D("pion_xbj_ID", "Reconstruction x PID efficiency  |  x_B  |  18x275 GeV", nbins, log_bins_xbj.data());
    pion_xbj_ID->Divide(pion23, pion21);
    pion_xbj_ID->SetLineColor(kRed);
    pion_xbj_ID->SetStats(kFALSE);
    pion_xbj_ID->GetXaxis()->SetTitle("x_B");
    pion_xbj_ID->GetYaxis()->SetTitle("Rec / MC events");
    pion_xbj_ID->GetYaxis()->SetRangeUser(0, 1.2);
    TH1D *kaon_xbj_ID = new TH1D("kaon_xbj_ID", "Reconstruction x PID efficiency  |  x_B  |  18x275 GeV", nbins, log_bins_xbj.data());
    kaon_xbj_ID->Divide(kaon23,kaon21);
    kaon_xbj_ID->SetLineColor(kBlue);
    kaon_xbj_ID->GetXaxis()->SetTitle("x_B");
    kaon_xbj_ID->GetYaxis()->SetTitle("Rec / MC events");
    kaon_xbj_ID->GetYaxis()->SetRangeUser(0, 1.2);
    kaon_xbj_ID->SetStats(kFALSE);

    TLegend *legend3 = new TLegend(0.17, 0.78, 0.30, 0.88);
    legend3->AddEntry(pion_xbj_ID, "Pions", "l");
    legend3->AddEntry(kaon_xbj_ID, "Kaons", "l");

    TCanvas *Ratio_xbj_ID = new TCanvas("Ratio_xbj_ID", "Ratio Canvas", 800, 600);
    Ratio_xbj_ID->SetLogx();
    kaon_xbj_ID->Draw("HIST");
    pion_xbj_ID->Draw("HIST SAME");
    //proton_xbj_ID->Draw("HIST SAME");
    legend3->Draw();
    
    Ratio_xbj_ID->Write();

    // FRACTION____________________________________________________________________________________________________
    TCanvas *MCfraction_xbj = new TCanvas("MCfraction_xbj", "Production with Q2", 800, 600);
    MCfraction_xbj->SetLogx();

    //int nBins = fpion21->GetNbinsX();
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = fpion21->GetBinContent(i) + fkaon21->GetBinContent(i) + fproton21->GetBinContent(i);

        if (sumBinContents > 0) {
            fpion21->SetBinContent(i, fpion21->GetBinContent(i) / (sumBinContents));
            fkaon21->SetBinContent(i, fkaon21->GetBinContent(i) / (sumBinContents));
            fproton21->SetBinContent(i, fproton21->GetBinContent(i) / (sumBinContents));
        }
    }
    fpion21->SetStats(kFALSE);
    fpion21->SetLineColor(kRed);
    fpion21->SetTitle("MC Production of charged particles | x_B | 18x275 GeV");
    fpion21->GetXaxis()->SetTitle("x_B");
    fpion21->GetYaxis()->SetRangeUser(0, 1.);
    fpion21->GetYaxis()->SetTitle("Relative fraction");
    fpion21->Draw("HIST");
    fkaon21->SetLineColor(kBlue);
    fkaon21->Draw("HIST SAME");
    fproton21->SetLineColor(kOrange);
    fproton21->Draw("HIST SAME");
    legend2->Draw();

    //MCfraction_xbj->SetTitle("Production of charged particles with Q2");
    MCfraction_xbj->Update();
    MCfraction_xbj->Write();

    // RECONSTRUCTED FRACTION
    TCanvas *RecFraction_xbj = new TCanvas("RecFraction_xbj", "Production with Q2", 800, 600);
    RecFraction_xbj->SetLogx();

    //int nBins = fpion22->GetNbinsX();
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = fpion22->GetBinContent(i) + fkaon22->GetBinContent(i) + fproton22->GetBinContent(i);

        if (sumBinContents > 0) {
            fpion22->SetBinContent(i, fpion22->GetBinContent(i) / (sumBinContents));
            fkaon22->SetBinContent(i, fkaon22->GetBinContent(i) / (sumBinContents));
            fproton22->SetBinContent(i, fproton22->GetBinContent(i) / (sumBinContents));
        }
    }
    fpion22->SetStats(kFALSE);
    fpion22->SetLineColor(kRed);
    fpion22->SetTitle("Reconstructed production of charged particles | x_B | 18x275 GeV");
    fpion22->GetXaxis()->SetTitle("x_B");
    fpion22->GetYaxis()->SetRangeUser(0, 1.);
    fpion22->GetYaxis()->SetTitle("Relative fraction");
    fpion22->Draw("HIST");
    fkaon22->SetLineColor(kBlue);
    fkaon22->Draw("HIST SAME");
    fproton22->SetLineColor(kOrange);
    fproton22->Draw("HIST SAME");
    legend2->Draw();

    //RecFraction_xbj->SetTitle("Production of charged particles with Q2");
    RecFraction_xbj->Update();
    RecFraction_xbj->Write();

    // ______________________________________________________________ z _____________________________________________________________

    double xmin_z = 1e-4;
    double xmax_z = 1;
    std::vector<double> log_bins_z = CreateLogBinning(nbins, xmin_z, xmax_z);
    TH1D *pion_z = new TH1D("pion_z", "Efficiency reconstruction with MC ID  |  z  |  18x275 GeV", nbins, log_bins_z.data());
    pion_z->Divide(pion32, pion31);
    pion_z->SetLineColor(kRed);
    pion_z->SetStats(kFALSE);
    pion_z->SetLineWidth(1);
    pion_z->GetXaxis()->SetTitle("z");
    pion_z->GetYaxis()->SetTitle("Rec / MC events");
    pion_z->GetYaxis()->SetRangeUser(0, 1.2);
    TH1D *kaon_z = new TH1D("kaon_z", "Ratio of Histograms", nbins, log_bins_z.data());
    kaon_z->Divide(kaon32,kaon31);
    kaon_z->SetLineColor(kBlue);
    TH1D *proton_z = new TH1D("proton_z", "Ratio of Histograms", nbins, log_bins_z.data());
    proton_z->Divide(proton32,proton31);
    proton_z->SetLineColor(kOrange);

    TCanvas *Ratio_z = new TCanvas("Ratio_z", "Ratio Canvas", 800, 600);
    Ratio_z->SetLogx();
    pion_z->Draw("HIST");
    kaon_z->Draw("HIST SAME");
    proton_z->Draw("HIST SAME");
    legend2->Draw();
    
    Ratio_z->Write();

    TH1D *pion_z_ID = new TH1D("pion_z_ID", "Reconstruction x PID efficiency  |  z  |  18x275 GeV", nbins, log_bins_z.data());
    pion_z_ID->Divide(pion33, pion31);
    pion_z_ID->SetLineColor(kRed);
    pion_z_ID->SetStats(kFALSE);
    //pion_z_ID->SetLineWidth(2);
    pion_z_ID->GetXaxis()->SetTitle("z");
    pion_z_ID->GetYaxis()->SetTitle("Rec / MC events");
    pion_z_ID->GetYaxis()->SetRangeUser(0, 1.2);
    TH1D *kaon_z_ID = new TH1D("kaon_z_ID", "Reconstruction x PID efficiency  |  z  |  18x275 GeV", nbins, log_bins_z.data());
    kaon_z_ID->Divide(kaon33,kaon31);
    kaon_z_ID->SetLineColor(kBlue);
    kaon_z_ID->SetStats(kFALSE);
    kaon_z_ID->GetXaxis()->SetTitle("z");
    kaon_z_ID->GetYaxis()->SetTitle("Rec / MC events");
    kaon_z_ID->GetYaxis()->SetRangeUser(0, 1.2);

    TCanvas *Ratio_z_ID = new TCanvas("Ratio_z_ID", "Ratio Canvas", 800, 600);
    Ratio_z_ID->SetLogx();
    kaon_z_ID->Draw("HIST");
    pion_z_ID->Draw("HIST SAME");
    //proton_z_ID->Draw("HIST SAME");
    legend->Draw();
    
    Ratio_z_ID->Write();

    // FRACTION____________________________________________________________________________________________________
    TCanvas *MCfraction_z = new TCanvas("MCfraction_z", "Production with Q2", 800, 600);
    MCfraction_z->SetLogx();

    //int nBins = fpion31->GetNbinsX();
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = fpion31->GetBinContent(i) + fkaon31->GetBinContent(i) + fproton31->GetBinContent(i);

        if (sumBinContents > 0) {
            fpion31->SetBinContent(i, fpion31->GetBinContent(i) / (sumBinContents));
            fkaon31->SetBinContent(i, fkaon31->GetBinContent(i) / (sumBinContents));
            fproton31->SetBinContent(i, fproton31->GetBinContent(i) / (sumBinContents));
        }
    }
    fpion31->SetStats(kFALSE);
    fpion31->SetLineColor(kRed);
    fpion31->SetTitle("MC Production of charged particles | z | 18x275 GeV");
    fpion31->GetXaxis()->SetTitle("z");
    fpion31->GetYaxis()->SetRangeUser(0, 1.);
    fpion31->GetYaxis()->SetTitle("Relative fraction");
    fpion31->Draw("HIST");
    fkaon31->SetLineColor(kBlue);
    fkaon31->Draw("HIST SAME");
    fproton31->SetLineColor(kOrange);
    fproton31->Draw("HIST SAME");
    legend2->Draw();

    //MCfraction_z->SetTitle("Production of charged particles with Q2");
    MCfraction_z->Update();
    MCfraction_z->Write();

    // RECONSTRUCTED FRACTION
    TCanvas *RecFraction_z = new TCanvas("RecFraction_z", "Production with Q2", 800, 600);
    RecFraction_z->SetLogx();

    //int nBins = fpion32->GetNbinsX();
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = fpion32->GetBinContent(i) + fkaon32->GetBinContent(i) + fproton32->GetBinContent(i);

        if (sumBinContents > 0) {
            fpion32->SetBinContent(i, fpion32->GetBinContent(i) / (sumBinContents));
            fkaon32->SetBinContent(i, fkaon32->GetBinContent(i) / (sumBinContents));
            fproton32->SetBinContent(i, fproton32->GetBinContent(i) / (sumBinContents));
        }
    }
    fpion32->SetStats(kFALSE);
    fpion32->SetLineColor(kRed);
    fpion32->SetTitle("Reconstructed production of charged particles | z | 18x275 GeV");
    fpion32->GetXaxis()->SetTitle("z");
    fpion32->GetYaxis()->SetRangeUser(0, 1.);
    fpion32->GetYaxis()->SetTitle("Relative fraction");
    fpion32->Draw("HIST");
    fkaon32->SetLineColor(kBlue);
    fkaon32->Draw("HIST SAME");
    fproton32->SetLineColor(kOrange);
    fproton32->Draw("HIST SAME");
    legend2->Draw();

    //RecFraction_z->SetTitle("Production of charged particles with Q2");
    RecFraction_z->Update();
    RecFraction_z->Write();


    // ______________________________________________________________ P_hT _____________________________________________________________

    double xmin_PhT= 1e-2;
    double xmax_PhT = 50;
    std::vector<double> log_bins_PhT = CreateLogBinning(nbins, xmin_PhT, xmax_PhT);
    TH1D *pion_PhT = new TH1D("pion_PhT", "Efficiency reconstruction with MC ID  |  P_hT  |  18x275 GeV", nbins, log_bins_PhT.data());
    pion_PhT->Divide(pion42, pion41);
    pion_PhT->SetLineColor(kRed);
    pion_PhT->SetStats(kFALSE);
    pion_PhT->SetLineWidth(1);
    pion_PhT->GetXaxis()->SetTitle("P_hT [GeV]");
    pion_PhT->GetYaxis()->SetTitle("Rec / MC events");
    pion_PhT->GetYaxis()->SetRangeUser(0, 1.2);
    TH1D *kaon_PhT = new TH1D("kaon_PhT", "Ratio of Histograms", nbins, log_bins_PhT.data());
    kaon_PhT->Divide(kaon42,kaon41);
    kaon_PhT->SetLineColor(kBlue);
    TH1D *proton_PhT = new TH1D("proton_PhT", "Ratio of Histograms", nbins, log_bins_PhT.data());
    proton_PhT->Divide(proton42,proton41);
    proton_PhT->SetLineColor(kOrange);

    TLegend *legend1 = new TLegend(0.17, 0.78, 0.30, 0.88);
    legend1->AddEntry(pion_PhT, "Pions", "l");
    legend1->AddEntry(kaon_PhT, "Kaons", "l");
    legend1->AddEntry(proton_PhT, "Protons", "l");

    TCanvas *Ratio_PhT = new TCanvas("Ratio_PhT", "Ratio Canvas", 800, 600);
    Ratio_PhT->SetLogx();
    pion_PhT->Draw("HIST");
    kaon_PhT->Draw("HIST SAME");
    proton_PhT->Draw("HIST SAME");
    legend1->Draw();
    
    Ratio_PhT->Write();

    TH1D *pion_PhT_ID = new TH1D("pion_PhT_ID", "Reconstruction x PID efficiency  |  P_hT  |  18x275 GeV", nbins, log_bins_PhT.data());
    pion_PhT_ID->Divide(pion43, pion41);
    pion_PhT_ID->SetLineColor(kRed);
    pion_PhT_ID->SetStats(kFALSE);
    //pion_PhT_ID->SetLineWidth(2);
    pion_PhT_ID->GetXaxis()->SetTitle("P_hT [GeV]");
    pion_PhT_ID->GetYaxis()->SetTitle("Rec / MC events");
    pion_PhT_ID->GetYaxis()->SetRangeUser(0, 1.2);
    TH1D *kaon_PhT_ID = new TH1D("kaon_PhT_ID", "Reconstruction x PID efficiency  |  P_hT  |  18x275 GeV", nbins, log_bins_PhT.data());
    kaon_PhT_ID->Divide(kaon43,kaon41);
    kaon_PhT_ID->GetXaxis()->SetTitle("P_hT [GeV]");
    kaon_PhT_ID->GetYaxis()->SetTitle("Rec / MC events");
    kaon_PhT_ID->SetLineColor(kBlue);
    kaon_PhT_ID->SetStats(kFALSE);
    kaon_PhT_ID->GetYaxis()->SetRangeUser(0, 1.2);

    TCanvas *Ratio_PhT_ID = new TCanvas("Ratio_PhT_ID", "Ratio Canvas", 800, 600);
    Ratio_PhT_ID->SetLogx();
    kaon_PhT_ID->Draw("HIST");
    pion_PhT_ID->Draw("HIST SAME");
    //proton_PhT_ID->Draw("HIST SAME");
    legend3->Draw();
    
    Ratio_PhT_ID->Write();

    // FRACTION____________________________________________________________________________________________________
    TCanvas *MCfraction_PhT = new TCanvas("MCfraction_PhT", "Production with Q2", 800, 600);
    MCfraction_PhT->SetLogx();

    //int nBins = fpion41->GetNbinsX();
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = fpion41->GetBinContent(i) + fkaon41->GetBinContent(i) + fproton41->GetBinContent(i);

        if (sumBinContents > 0) {
            fpion41->SetBinContent(i, fpion41->GetBinContent(i) / (sumBinContents));
            fkaon41->SetBinContent(i, fkaon41->GetBinContent(i) / (sumBinContents));
            fproton41->SetBinContent(i, fproton41->GetBinContent(i) / (sumBinContents));
        }
    }
    fpion41->SetStats(kFALSE);
    fpion41->SetLineColor(kRed);
    fpion41->SetTitle("MC Production of charged particles | P_hT | 18x275 GeV");
    fpion41->GetXaxis()->SetTitle("P_hT [GeV]");
    fpion41->GetYaxis()->SetRangeUser(0, 1.);
    fpion41->GetYaxis()->SetTitle("Relative fraction");
    fpion41->Draw("HIST");
    fkaon41->SetLineColor(kBlue);
    fkaon41->Draw("HIST SAME");
    fproton41->SetLineColor(kOrange);
    fproton41->Draw("HIST SAME");
    legend2->Draw();

    //MCfraction_PhT->SetTitle("Production of charged particles with Q2");
    MCfraction_PhT->Update();
    MCfraction_PhT->Write();

    // RECONSTRUCTED FRACTION
    TCanvas *RecFraction_PhT = new TCanvas("RecFraction_PhT", "Production with Q2", 800, 600);
    RecFraction_PhT->SetLogx();

    //int nBins = fpion42->GetNbinsX();
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBins; ++i) {
        double sumBinContents = fpion42->GetBinContent(i) + fkaon42->GetBinContent(i) + fproton42->GetBinContent(i);

        if (sumBinContents > 0) {
            fpion42->SetBinContent(i, fpion42->GetBinContent(i) / (sumBinContents));
            fkaon42->SetBinContent(i, fkaon42->GetBinContent(i) / (sumBinContents));
            fproton42->SetBinContent(i, fproton42->GetBinContent(i) / (sumBinContents));
        }
    }
    fpion42->SetStats(kFALSE);
    fpion42->SetLineColor(kRed);
    fpion42->SetTitle("Reconstructed production of charged particles | P_hT | 18x275 GeV");
    fpion42->GetXaxis()->SetTitle("P_hT [GeV]");
    fpion42->GetYaxis()->SetRangeUser(0, 1.);
    fpion42->GetYaxis()->SetTitle("Relative fraction");
    fpion42->Draw("HIST");
    fkaon42->SetLineColor(kBlue);
    fkaon42->Draw("HIST SAME");
    fproton42->SetLineColor(kOrange);
    fproton42->Draw("HIST SAME");
    legend2->Draw();

    //RecFraction_PhT->SetTitle("Production of charged particles with Q2");
    RecFraction_PhT->Update();
    RecFraction_PhT->Write();



    // ______________________________________________________________ Eta _____________________________________________________________

    TH1D *pionEta = new TH1D("pionEta", "Efficiency reconstruction with MC ID  |  Eta  |  18x275 GeV", 80, 1.5, 3.5);
    pionEta->Divide(pion52, pion51);
    pionEta->SetLineColor(kRed);
    pionEta->SetStats(kFALSE);
    pionEta->GetXaxis()->SetTitle("Eta");
    pionEta->GetYaxis()->SetTitle("Rec / MC events");
    pionEta->GetYaxis()->SetRangeUser(0, 1.2);
    TH1D *kaonEta = new TH1D("kaonEta", "Ratio of Histograms", 80, 1.5, 3.5);
    kaonEta->Divide(kaon52,kaon51);
    kaonEta->SetLineColor(kBlue);
    TH1D *protonEta = new TH1D("protonEta", "Ratio of Histograms", 80, 1.5, 3.5);
    protonEta->Divide(proton52,proton51);
    protonEta->SetLineColor(kOrange);

    TCanvas *Ratio_Eta = new TCanvas("Ratio_Eta", "Ratio Canvas", 800, 600);
    pionEta->Draw("HIST");
    kaonEta->Draw("HIST SAME");
    protonEta->Draw("HIST SAME");
    legend1->Draw();
    
    Ratio_Eta->Write();

    TH1D *pionEta_ID = new TH1D("pionEta_ID", "Reconstruction x PID efficiency  |  Eta  |  18x275 GeV", 80, 1.5, 3.5);
    pionEta_ID->Divide(pion53, pion51);
    pionEta_ID->SetLineColor(kRed);
    pionEta_ID->SetStats(kFALSE);
    //pionEta_ID->SetLineWidth(2);
    pionEta_ID->GetXaxis()->SetTitle("Eta");
    pionEta_ID->GetYaxis()->SetTitle("Rec / MC events");
    pionEta_ID->GetYaxis()->SetRangeUser(0, 1.2);
    TH1D *kaonEta_ID = new TH1D("kaonEta_ID", "Reconstruction x PID efficiency  |  Eta  |  18x275 GeV", 80, 1.5, 3.5);
    kaonEta_ID->Divide(kaon53,kaon51);
    kaonEta_ID->SetLineColor(kBlue);
    kaonEta_ID->SetStats(kFALSE);
    kaonEta_ID->GetXaxis()->SetTitle("Eta");
    kaonEta_ID->GetYaxis()->SetTitle("Rec / MC events");
    kaonEta_ID->GetYaxis()->SetRangeUser(0, 1.2);
    
    TCanvas *Ratio_Eta_ID = new TCanvas("Ratio_Eta_ID", "Ratio Canvas", 800, 600);
    kaonEta_ID->Draw("HIST");
    pionEta_ID->Draw("HIST SAME");
    //protonEta_ID->Draw("HIST SAME");
    legend3->Draw();
    
    Ratio_Eta_ID->Write();

    // FRACTION____________________________________________________________________________________________________
    TCanvas *MCfraction_Eta = new TCanvas("MCfraction_Eta", "Production with Q2", 800, 600);
    //MCfraction_Eta->SetLogx();

    int nBin = fpion51->GetNbinsX();
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBin; ++i) {
        double sumBinContents = fpion51->GetBinContent(i) + fkaon51->GetBinContent(i) + fproton51->GetBinContent(i);

        if (sumBinContents > 0) {
            fpion51->SetBinContent(i, fpion51->GetBinContent(i) / (sumBinContents));
            fkaon51->SetBinContent(i, fkaon51->GetBinContent(i) / (sumBinContents));
            fproton51->SetBinContent(i, fproton51->GetBinContent(i) / (sumBinContents));
        }
    }
    fpion51->SetStats(kFALSE);
    fpion51->SetLineColor(kRed);
    fpion51->SetTitle("MC Production of charged particles | Eta | 18x275 GeV");
    fpion51->GetXaxis()->SetTitle("Eta");
    fpion51->GetYaxis()->SetRangeUser(0, 1.);
    fpion51->GetYaxis()->SetTitle("Relative fraction");
    fpion51->Draw("HIST");
    fkaon51->SetLineColor(kBlue);
    fkaon51->Draw("HIST SAME");
    fproton51->SetLineColor(kOrange);
    fproton51->Draw("HIST SAME");
    legend2->Draw();

    //MCfraction_Eta->SetTitle("Production of charged particles with Q2");
    MCfraction_Eta->Update();
    MCfraction_Eta->Write();

    // RECONSTRUCTED FRACTION
    TCanvas *RecFraction_Eta = new TCanvas("RecFraction_Eta", "Production with Q2", 800, 600);
    //RecFraction_Eta->SetLogx();

    //int nBins = fpion52->GetNbinsX();
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBin; ++i) {
        double sumBinContents = fpion52->GetBinContent(i) + fkaon52->GetBinContent(i) + fproton52->GetBinContent(i);

        if (sumBinContents > 0) {
            fpion52->SetBinContent(i, fpion52->GetBinContent(i) / (sumBinContents));
            fkaon52->SetBinContent(i, fkaon52->GetBinContent(i) / (sumBinContents));
            fproton52->SetBinContent(i, fproton52->GetBinContent(i) / (sumBinContents));
        }
    }
    fpion52->SetStats(kFALSE);
    fpion52->SetLineColor(kRed);
    fpion52->SetTitle("Reconstructed production of charged particles | Eta | 18x275 GeV");
    fpion52->GetXaxis()->SetTitle("Eta");
    fpion52->GetYaxis()->SetRangeUser(0, 1.);
    fpion52->GetYaxis()->SetTitle("Relative fraction");
    fpion52->Draw("HIST");
    fkaon52->SetLineColor(kBlue);
    fkaon52->Draw("HIST SAME");
    fproton52->SetLineColor(kOrange);
    fproton52->Draw("HIST SAME");
    legend2->Draw();

    //RecFraction_Eta->SetTitle("Production of charged particles with Q2");
    RecFraction_Eta->Update();
    RecFraction_Eta->Write();


    // ______________________________________________________________ Phi _____________________________________________________________

    TH1D *pionPhi = new TH1D("pionPhi", "Efficiency reconstruction with MC ID  |  Phi  |  18x275 GeV", 80, 0, 25);
    pionPhi->Divide(pion62, pion61);
    pionPhi->SetLineColor(kRed);
    pionPhi->SetStats(kFALSE);
    pionPhi->GetXaxis()->SetTitle("Phi [Deg]");
    pionPhi->GetYaxis()->SetTitle("Rec / MC events");
    pionPhi->GetYaxis()->SetRangeUser(0, 1.2);
    TH1D *kaonPhi = new TH1D("kaonPhi", "Ratio of Histograms", 80,  0, 25);
    kaonPhi->Divide(kaon62,kaon61);
    kaonPhi->SetLineColor(kBlue);
    TH1D *protonPhi = new TH1D("protonPhi", "Ratio of Histograms", 80,  0, 25);
    protonPhi->Divide(proton62,proton61);
    protonPhi->SetLineColor(kOrange);

    TCanvas *Ratio_Phi = new TCanvas("Ratio_Phi", "Ratio Canvas", 800, 600);
    pionPhi->Draw("HIST");
    kaonPhi->Draw("HIST SAME");
    protonPhi->Draw("HIST SAME");
    legend2->Draw();
    
    Ratio_Phi->Write();

    TH1D *pionPhi_ID = new TH1D("pionPhi_ID", "Reconstruction x PID efficiency  |  Phi  |  18x275 GeV", 80, 0, 25);
    pionPhi_ID->Divide(pion63, pion61);
    pionPhi_ID->SetLineColor(kRed);
    pionPhi_ID->SetStats(kFALSE);
    pionPhi_ID->GetXaxis()->SetTitle("Phi [Deg]");
    pionPhi_ID->GetYaxis()->SetTitle("Rec / MC events");
    pionPhi_ID->GetYaxis()->SetRangeUser(0, 1.2);
    TH1D *kaonPhi_ID = new TH1D("kaonPhi_ID", "Reconstruction x PID efficiency  |  Phi  |  18x275 GeV", 80,  0, 25);
    kaonPhi_ID->Divide(kaon63,kaon61);
    kaonPhi_ID->SetLineColor(kBlue);
    kaonPhi_ID->SetStats(kFALSE);
    kaonPhi_ID->GetXaxis()->SetTitle("Phi [Deg]");
    kaonPhi_ID->GetYaxis()->SetTitle("Rec / MC events");
    kaonPhi_ID->GetYaxis()->SetRangeUser(0, 1.2);

    TCanvas *Ratio_Phi_ID = new TCanvas("Ratio_Phi_ID", "Ratio Canvas", 800, 600);
    kaonPhi_ID->Draw("HIST");
    pionPhi_ID->Draw("HIST SAME");
    legend->Draw();
    
    Ratio_Phi_ID->Write();

    // FRACTION____________________________________________________________________________________________________
    TCanvas *MCfraction_Phi = new TCanvas("MCfraction_Phi", "Production with Q2", 800, 600);
    //MCfraction_Phi->SetLogx();

    //int nBin = fpion61->GetNbinsX();
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBin; ++i) {
        double sumBinContents = fpion61->GetBinContent(i) + fkaon61->GetBinContent(i) + fproton61->GetBinContent(i);

        if (sumBinContents > 0) {
            fpion61->SetBinContent(i, fpion61->GetBinContent(i) / (sumBinContents));
            fkaon61->SetBinContent(i, fkaon61->GetBinContent(i) / (sumBinContents));
            fproton61->SetBinContent(i, fproton61->GetBinContent(i) / (sumBinContents));
        }
    }
    fpion61->SetStats(kFALSE);
    fpion61->SetLineColor(kRed);
    fpion61->SetTitle("MC Production of charged particles | Phi | 18x275 GeV");
    fpion61->GetXaxis()->SetTitle("Phi");
    fpion61->GetYaxis()->SetRangeUser(0, 1.);
    fpion61->GetYaxis()->SetTitle("Relative fraction");
    fpion61->Draw("HIST");
    fkaon61->SetLineColor(kBlue);
    fkaon61->Draw("HIST SAME");
    fproton61->SetLineColor(kOrange);
    fproton61->Draw("HIST SAME");
    legend2->Draw();

    //MCfraction_Phi->SetTitle("Production of charged particles with Q2");
    MCfraction_Phi->Update();
    MCfraction_Phi->Write();

    // RECONSTRUCTED FRACTION
    TCanvas *RecFraction_Phi = new TCanvas("RecFraction_Phi", "Production with Q2", 800, 600);
    //RecFraction_Phi->SetLogx();

    //int nBins = fpion62->GetNbinsX();
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBin; ++i) {
        double sumBinContents = fpion62->GetBinContent(i) + fkaon62->GetBinContent(i) + fproton62->GetBinContent(i);

        if (sumBinContents > 0) {
            fpion62->SetBinContent(i, fpion62->GetBinContent(i) / (sumBinContents));
            fkaon62->SetBinContent(i, fkaon62->GetBinContent(i) / (sumBinContents));
            fproton62->SetBinContent(i, fproton62->GetBinContent(i) / (sumBinContents));
        }
    }
    fpion62->SetStats(kFALSE);
    fpion62->SetLineColor(kRed);
    fpion62->SetTitle("Reconstructed production of charged particles | Phi | 18x275 GeV");
    fpion62->GetXaxis()->SetTitle("Phi");
    fpion62->GetYaxis()->SetRangeUser(0, 1.);
    fpion62->GetYaxis()->SetTitle("Relative fraction");
    fpion62->Draw("HIST");
    fkaon62->SetLineColor(kBlue);
    fkaon62->Draw("HIST SAME");
    fproton62->SetLineColor(kOrange);
    fproton62->Draw("HIST SAME");
    legend2->Draw();

    //RecFraction_Phi->SetTitle("Production of charged particles with Q2");
    RecFraction_Phi->Update();
    RecFraction_Phi->Write();


    // ______________________________________________________________ Mom _____________________________________________________________

    double xminMom = 1e-1;
    double xmaxMom = 50;
    std::vector<double> log_binsMom = CreateLogBinning(nbins, xminMom, xmaxMom);
    TH1D *pionMom = new TH1D("pionMom", "Efficiency reconstruction with MC ID  |  P_h  |  18x275 GeV", nbins, log_binsMom.data());
    pionMom->Divide(pion72, pion71);
    pionMom->SetLineColor(kRed);
    pionMom->SetStats(kFALSE);
    pionMom->SetLineWidth(1);
    pionMom->GetXaxis()->SetTitle("P_h [GeV]");
    pionMom->GetYaxis()->SetTitle("Rec / MC events");
    pionMom->GetYaxis()->SetRangeUser(0, 1.2);
    TH1D *kaonMom = new TH1D("kaonMom", "Ratio of Histograms", nbins, log_binsMom.data());
    kaonMom->Divide(kaon72,kaon71);
    kaonMom->SetLineColor(kBlue);
    TH1D *protonMom = new TH1D("protonMom", "Ratio of Histograms", nbins, log_binsMom.data());
    protonMom->Divide(proton72,proton71);
    protonMom->SetLineColor(kOrange);

    TCanvas *Ratio_Mom = new TCanvas("Ratio_Mom", "Ratio Canvas", 800, 600);
    Ratio_Mom->SetLogx();
    pionMom->Draw("HIST");
    kaonMom->Draw("HIST SAME");
    protonMom->Draw("HIST SAME");
    legend1->Draw();
    
    Ratio_Mom->Write();

    TH1D *pionMom_ID = new TH1D("pionMom_ID", "Reconstruction x PID efficiency  |  P_h  |  18x275 GeV", nbins, log_binsMom.data());
    pionMom_ID->Divide(pion73, pion71);
    pionMom_ID->SetLineColor(kRed);
    pionMom_ID->SetStats(kFALSE);
    //pionMom_ID->SetLineWidth(2);
    pionMom_ID->GetXaxis()->SetTitle("P_h [GeV]");
    pionMom_ID->GetYaxis()->SetTitle("Rec / MC events");
    pionMom_ID->GetYaxis()->SetRangeUser(0, 1.2);
    TH1D *kaonMom_ID = new TH1D("kaonMom_ID", "Reconstruction x PID efficiency  |  P_h  |  18x275 GeV", nbins, log_binsMom.data());
    kaonMom_ID->Divide(kaon73,kaon71);
    kaonMom_ID->SetLineColor(kBlue);
    kaonMom_ID->SetStats(kFALSE);
    kaonMom_ID->GetXaxis()->SetTitle("P_h [GeV]");
    kaonMom_ID->GetYaxis()->SetTitle("Rec / MC events");
    kaonMom_ID->GetYaxis()->SetRangeUser(0, 1.2);

    TCanvas *Ratio_Mom_ID = new TCanvas("Ratio_Mom_ID", "Ratio Canvas", 800, 600);
    Ratio_Mom_ID->SetLogx();
    kaonMom_ID->Draw("HIST");
    pionMom_ID->Draw("HIST SAME");
    //protonMom_ID->Draw("HIST SAME");
    legend3->Draw();
    
    Ratio_Mom_ID->Write();

    // FRACTION____________________________________________________________________________________________________
    TCanvas *MCfraction_Mom = new TCanvas("MCfraction_Mom", "Production with Q2", 800, 600);
    MCfraction_Mom->SetLogx();

    //int nBin = fpion71->GetNbinsX();
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBin; ++i) {
        double sumBinContents = fpion71->GetBinContent(i) + fkaon71->GetBinContent(i) + fproton71->GetBinContent(i);

        if (sumBinContents > 0) {
            fpion71->SetBinContent(i, fpion71->GetBinContent(i) / (sumBinContents));
            fkaon71->SetBinContent(i, fkaon71->GetBinContent(i) / (sumBinContents));
            fproton71->SetBinContent(i, fproton71->GetBinContent(i) / (sumBinContents));
        }
    }
    fpion71->SetStats(kFALSE);
    fpion71->SetLineColor(kRed);
    fpion71->SetTitle("MC Production of charged particles | P_h | 18x275 GeV");
    fpion71->GetXaxis()->SetTitle("P_h [GeV]");
    fpion71->GetYaxis()->SetRangeUser(0, 1.);
    fpion71->GetYaxis()->SetTitle("Relative fraction");
    fpion71->Draw("HIST");
    fkaon71->SetLineColor(kBlue);
    fkaon71->Draw("HIST SAME");
    fproton71->SetLineColor(kOrange);
    fproton71->Draw("HIST SAME");
    legend2->Draw();

    //MCfraction_Mom->SetTitle("Production of charged particles with Q2");
    MCfraction_Mom->Update();
    MCfraction_Mom->Write();

    // RECONSTRUCTED FRACTION
    TCanvas *RecFraction_Mom = new TCanvas("RecFraction_Mom", "Production with Q2", 800, 600);
    RecFraction_Mom->SetLogx();

    //int nBins = fpion72->GetNbinsX();
    // Normalize each bin across all histograms
    for (int i = 1; i <= nBin; ++i) {
        double sumBinContents = fpion72->GetBinContent(i) + fkaon72->GetBinContent(i) + fproton72->GetBinContent(i);

        if (sumBinContents > 0) {
            fpion72->SetBinContent(i, fpion72->GetBinContent(i) / (sumBinContents));
            fkaon72->SetBinContent(i, fkaon72->GetBinContent(i) / (sumBinContents));
            fproton72->SetBinContent(i, fproton72->GetBinContent(i) / (sumBinContents));
        }
    }
    fpion72->SetStats(kFALSE);
    fpion72->SetLineColor(kRed);
    fpion72->SetTitle("Reconstructed production of charged particles | P_h | 18x275 GeV");
    fpion72->GetXaxis()->SetTitle("P_h [GeV]");
    fpion72->GetYaxis()->SetRangeUser(0, 1.);
    fpion72->GetYaxis()->SetTitle("Relative fraction");
    fpion72->Draw("HIST");
    fkaon72->SetLineColor(kBlue);
    fkaon72->Draw("HIST SAME");
    fproton72->SetLineColor(kOrange);
    fproton72->Draw("HIST SAME");
    legend2->Draw();

    //RecFraction_Mom->SetTitle("Production of charged particles with Q2");
    RecFraction_Mom->Update();
    RecFraction_Mom->Write();

    outputFile->Close();
    combinedFile->Close();

    return 0;
}
