#include <iostream>
#include <cmath>
#include <TStyle.h>
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"

void FourierAnalysis(const char* inputFileName, const char* outputFileName) {
    // Apertura del file ROOT di input
    TFile* inputFile = TFile::Open(inputFileName);
    //if (inputFile) inputFile->ls();
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Errore: impossibile aprire il file " << inputFileName << std::endl;
        return;
    }

    TCanvas* canvas = (TCanvas*)inputFile->Get("Fourier");
    if (!canvas) {
        std::cerr << "Errore: canvas 'Fourier' non trovata nel file " << inputFileName << std::endl;
        inputFile->Close();
        return;
    }

    TCanvas* canvas2 = (TCanvas*)inputFile->Get("Fourier2");
    if (!canvas2) {
        std::cerr << "Errore: canvas 'Fourier2' non trovata nel file " << inputFileName << std::endl;
        inputFile->Close();
        return;
    }

    TH1D* h_phi_sin = (TH1D*)inputFile->Get("DeltaPhi_Sin"); // Usa il nome esatto dell'istogramma
    if (!h_phi_sin) {
        std::cerr << "Errore: impossibile trovare l'istogramma 'DeltaPhi_Sin'." << std::endl;
        inputFile->Close();
        return;
    }

    TH1D* h_phi_xlow = (TH1D*)inputFile->Get("DeltaPhi_xLow"); // Usa il nome esatto dell'istogramma
    if (!h_phi_xlow) {
        std::cerr << "Errore: impossibile trovare l'istogramma 'DeltaPhi_xLow'." << std::endl;
        inputFile->Close();
        return;
    }

    TH1D* h_phi_xhigh = (TH1D*)inputFile->Get("DeltaPhi_xHigh"); // Usa il nome esatto dell'istogramma
    if (!h_phi_xhigh) {
        std::cerr << "Errore: impossibile trovare l'istogramma 'DeltaPhi_xHigh'." << std::endl;
        inputFile->Close();
        return;
    }

    // Recupera l'istogramma dalla canvas
    TH1F* h_phi = (TH1F*)canvas->GetPrimitive("DeltaPhi");
    //canvas->ls();
    if (!h_phi) {
        std::cerr << "Errore: impossibile trovare l'istogramma nella canvas 'Fourier'." << std::endl;
        inputFile->Close();
        return;
    }

    if (h_phi->GetEntries() == 0) {
        std::cerr << "Errore: l'istogramma uno è vuoto." << std::endl;
        inputFile->Close();
        return;
    }

    TH1F* h_phi2 = (TH1F*)canvas2->GetPrimitive("DeltaPhi2");
    //canvas->ls();
    if (!h_phi2) {
        std::cerr << "Errore: impossibile trovare l'istogramma nella canvas 'Fourier2'." << std::endl;
        inputFile->Close();
        return;
    }

    if (h_phi2->GetEntries() == 0) {
        std::cerr << "Errore: l'istogramma due è vuoto." << std::endl;
        inputFile->Close();
        return;
    }

    if (h_phi_sin->GetEntries() == 0) {
    std::cerr << "Errore: l'istogramma 'DeltaPhi_Sin' è vuoto." << std::endl;
    inputFile->Close();
    return;
    }

    // NORMALIZZO L'ISTOGRAMMA
    h_phi->Scale(10.0 / h_phi->Integral());
    h_phi2->Scale(1.0 / h_phi2->Integral());
    h_phi_sin->Scale(1.0 / h_phi_sin->Integral());
    h_phi_xlow->Scale(1.0 / h_phi_xlow->Integral());
    h_phi_xhigh->Scale(1.0 / h_phi_xhigh->Integral());

    // Creazione del file ROOT di output
    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Errore: impossibile creare il file " << outputFileName << std::endl;
        inputFile->Close();
        return;
    }

    TPaveStats* stats = (TPaveStats*)h_phi->FindObject("stats");
    if (stats) {
        stats->Delete(); // Rimuove il box delle statistiche precedente
    }

    TPaveStats* stats2 = (TPaveStats*)h_phi2->FindObject("stats");
    if (stats2) {
        stats2->Delete(); // Rimuove il box delle statistiche precedente
    }

    //__________________________________________________ FOURIER _____________________________________________________________________

    double a0 = 0.0, a1 = 0.0, b1 = 0.0, a2 = 0.0, b2 = 0.0, a3 = 0.0, b3 = 0.0;
    double deltaPhi = h_phi->GetBinWidth(1); // Larghezza del bin
    double normFactor = deltaPhi / TMath::Pi(); // Normalizzazione (per 2pi)

    int nBins = h_phi->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
        double phi = h_phi->GetBinCenter(i);
        double f_phi = h_phi->GetBinContent(i);

        a0 += f_phi * normFactor;
        a1 += f_phi * TMath::Cos(phi) * normFactor;
        b1 += f_phi * TMath::Sin(phi) * normFactor;
        a2 += f_phi * TMath::Cos(2*phi) * normFactor;
        b2 += f_phi * TMath::Sin(2*phi) * normFactor;
        a3 += f_phi * TMath::Cos(3*phi) * normFactor;
        b3 += f_phi * TMath::Sin(3*phi) * normFactor;
    }

    double a00 = 0.0, a11 = 0.0, b11 = 0.0, a22 = 0.0, b22 = 0.0, a33 = 0.0, b33 = 0.0;
    double deltaPhi2 = h_phi2->GetBinWidth(1); // Larghezza del bin
    double normFactor2 = deltaPhi2 / TMath::Pi(); // Normalizzazione (per 2pi)

    int nBins2 = h_phi2->GetNbinsX();
    for (int j = 1; j <= nBins2; ++j) {
        double phi2 = h_phi2->GetBinCenter(j);
        double f_phi2 = h_phi2->GetBinContent(j);

        a00 += f_phi2 * normFactor2;
        a11 += f_phi2 * TMath::Cos(phi2) * normFactor2;
        b11 += f_phi2 * TMath::Sin(phi2) * normFactor2;
        a22 += f_phi2 * TMath::Cos(2*phi2) * normFactor2;
        b22 += f_phi2 * TMath::Sin(2*phi2) * normFactor2;
        a33 += f_phi2 * TMath::Cos(3*phi2) * normFactor2;
        b33 += f_phi2 * TMath::Sin(3*phi2) * normFactor2;
    }

    // Creazione del canvas
    TCanvas* c = new TCanvas("DeltaPhi", "Analisi di Fourier", 800, 600);
    c->ToggleEditor();
    // Creazione del fit
    TF1* fitFunc = new TF1("fitFunc", "[0] + [1]*cos(x) + [2]*sin(x) + [3]*cos(2*x) + [4]*sin(2*x) ", h_phi->GetXaxis()->GetXmin(), h_phi->GetXaxis()->GetXmax());
    fitFunc->SetParameters(a0, a1, b1, a2, b2);
    fitFunc->SetParName(0, "a0");
    fitFunc->SetParName(1, "a1");
    fitFunc->SetParName(2, "b1");
    fitFunc->SetParName(3, "a2");
    fitFunc->SetParName(4, "b2");
    fitFunc->SetLineColor(kRed);
    h_phi->Draw("E");
    h_phi->Fit(fitFunc, "R");
    fitFunc->Draw("SAME");
    gStyle->SetOptFit(111);
    TPaveStats* stats3 = (TPaveStats*)h_phi->GetListOfFunctions()->FindObject("stats");
    if (stats3) {
        stats3->SetX1NDC(0.12);     // Sposta il pannello
        stats3->SetX2NDC(0.38);     
        stats3->SetY1NDC(0.55);     
        stats3->SetY2NDC(0.88);     
    }
    gPad->Update();
    c->Update();
    c->Write(); // Scrive il canvas nel file di output

    TCanvas* c2 = new TCanvas("DeltaPhi2", "Analisi di Fourier", 800, 600);
    c2->ToggleEditor();
    // Creazione del fit
    TF1* fitFunc2 = new TF1("fitFunc2", "[0] + [1]*cos(x) + [2]*sin(x) + [3]*cos(2*x) + [4]*sin(2*x)", h_phi2->GetXaxis()->GetXmin(), h_phi2->GetXaxis()->GetXmax());
    fitFunc2->SetParameters(a00, a11, b11, a22, b22);
    fitFunc2->SetParName(0, "a0");
    fitFunc2->SetParName(1, "a1");
    fitFunc2->SetParName(2, "b1");
    fitFunc2->SetParName(3, "a2");
    fitFunc2->SetParName(4, "b2");
    fitFunc2->SetLineColor(kRed);
    h_phi2->Draw("E");
    h_phi2->Fit(fitFunc2, "R");
    fitFunc2->Draw("SAME");
    gStyle->SetOptFit(111);
    TPaveStats* stats4 = (TPaveStats*)h_phi2->GetListOfFunctions()->FindObject("stats");
    if (stats4) {
        stats4->SetX1NDC(0.12);     // Sposta il pannello
        stats4->SetX2NDC(0.38);     
        stats4->SetY1NDC(0.50);     
        stats4->SetY2NDC(0.88);     
    }
    gPad->Update();
    c2->Update();
    c2->Write();

    //__________________________________________________________ SIVERS ________________________________________________________________

    TCanvas* siv = new TCanvas("Sivers", "Sivers analysis", 800, 600);
    siv->ToggleEditor();
    TF1* fitFuncS = new TF1("fitFuncS", "[0]*sin(x) + [1] ", -TMath::Pi(), TMath::Pi());
    fitFuncS->SetParameter(0, 0.1);
    fitFuncS->SetParLimits(0, -1.0, 1.0); // Limita A_UT nell'intervallo fisico [-1, 1]
    fitFuncS->SetParName(0, "A_UT");
    fitFuncS->SetParName(1, "C");
    h_phi_sin->Draw("E");
    h_phi_sin->Fit(fitFuncS, "R");
    double A_UT = fitFuncS->GetParameter(0);
    double A_UT_err = fitFuncS->GetParError(0); 
    gStyle->SetOptFit(111);
    TPaveStats* statsS2 = (TPaveStats*)h_phi_sin->GetListOfFunctions()->FindObject("stats");
    if (statsS2) {
        statsS2->SetX1NDC(0.12);     // Sposta il pannello
        statsS2->SetX2NDC(0.35);     
        statsS2->SetY1NDC(0.62);     
        statsS2->SetY2NDC(0.88);     
    }
    gPad->Update();
    siv->Update();
    siv->Write();

    TCanvas* siv2 = new TCanvas("Sivers_xLow", "Sivers analysis", 800, 600);
    siv2->ToggleEditor();
    TF1* fitFuncS2 = new TF1("fitFuncS2", "[0]*sin(x) + [1]", -TMath::Pi(), TMath::Pi());
    fitFuncS2->SetParameter(0, 0.1);
    fitFuncS2->SetParLimits(0, -1.0, 1.0); // Limita A_UT nell'intervallo fisico [-1, 1]
    fitFuncS2->SetParName(0, "A_UT");
    fitFuncS2->SetParName(1, "C");
    h_phi_xlow->Draw("E");
    h_phi_xlow->Fit(fitFuncS2, "R");
    double A_UT_xlow = fitFuncS2->GetParameter(0);
    double A_UT_err_xlow = fitFuncS2->GetParError(0); 
    gStyle->SetOptFit(111);
    TPaveStats* statsS3 = (TPaveStats*)h_phi_xlow->GetListOfFunctions()->FindObject("stats");
    if (statsS3) {
        statsS3->SetX1NDC(0.12);     // Sposta il pannello
        statsS3->SetX2NDC(0.35);     
        statsS3->SetY1NDC(0.62);     
        statsS3->SetY2NDC(0.88);     
    }
    gPad->Update();
    siv2->Update();
    siv2->Write();

    TCanvas* siv3 = new TCanvas("Sivers_xHigh", "Sivers analysis", 800, 600);
    siv->ToggleEditor();
    TF1* fitFuncS3 = new TF1("fitFuncS3", "[0]*sin(x) + [1]", -TMath::Pi(), TMath::Pi());
    fitFuncS3->SetParameter(0, 0.1);
    fitFuncS3->SetParLimits(0, -1.0, 1.0); // Limita A_UT nell'intervallo fisico [-1, 1]
    fitFuncS3->SetParName(0, "A_UT");
    fitFuncS3->SetParName(1, "C");
    h_phi_xhigh->Draw("E");
    h_phi_xhigh->Fit(fitFuncS3, "R");
    double A_UT_xhigh = fitFuncS3->GetParameter(0);
    double A_UT_err_xhigh = fitFuncS3->GetParError(0); 
    gStyle->SetOptFit(111);
    TPaveStats* statsS4 = (TPaveStats*)h_phi_xhigh->GetListOfFunctions()->FindObject("stats");
    if (statsS4) {
        statsS4->SetX1NDC(0.12);     // Sposta il pannello
        statsS4->SetX2NDC(0.35);     
        statsS4->SetY1NDC(0.62);     
        statsS4->SetY2NDC(0.88);     
    }
    gPad->Update();
    siv3->Update();
    siv3->Write();

    // _________________________________________________________________________________________________________________________________
    
    // Pulizia
    delete c;
    delete c2;
    delete siv;
    delete siv2;
    delete siv3;
    outputFile->Write();
    outputFile->Close();
    inputFile->Close();

    std::cout << "Analisi completata. Risultati salvati in " << outputFileName << std::endl;
}

int main() {
    const char* inputFileName = "A11.combinedhistAbs(dRICHEta).root";
    const char* outputFileName = "fourierCoeff.root";
    FourierAnalysis(inputFileName, outputFileName);
    return 0;
}
