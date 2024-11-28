#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TKey.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TList.h>
#include <TObject.h>

void mergeHisto(const std::vector<std::string>& inputFiles, const char* outputFileName) {
    std::vector<TFile*> files;

    // Apri tutti i file di input
    for (const auto& inputFile : inputFiles) {
        TFile* file = TFile::Open(inputFile.c_str());
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Cannot open input file " << inputFile << std::endl;
            for (auto* f : files) {
                if (f) f->Close();
            }
            return;
        }
        files.push_back(file);
    }

    // Crea il file di output
    TFile* outputFile = TFile::Open(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << outputFileName << std::endl;
        for (auto* f : files) {
            if (f) f->Close();
        }
        return;
    }

    // Unisci i grafici e le canvas
    TList* keyList = files[0]->GetListOfKeys();
    TIter next(keyList);
    TKey* key;
    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        if (obj->InheritsFrom("TH1")) {
            TH1* hist = (TH1*)obj->Clone();
            hist->SetDirectory(nullptr);
            for (size_t i = 1; i < files.size(); ++i) {
                TH1* otherHist = (TH1*)files[i]->Get(hist->GetName());
                if (otherHist) {
                    hist->Add(otherHist);
                } else {
                    std::cerr << "Warning: Histogram " << hist->GetName()
                              << " not found in file " << inputFiles[i] << std::endl;
                }
            }
            outputFile->cd();
            hist->Write();
            delete hist;
        } else if (obj->InheritsFrom("TCanvas")) {
            // Clona la canvas di riferimento
            TCanvas* canvas = (TCanvas*)obj->Clone();
            canvas->SetName(Form("%s", canvas->GetName()));

            // Accedi ai primitivi e combinali
            for (size_t i = 1; i < files.size(); ++i) {
                TCanvas* otherCanvas = (TCanvas*)files[i]->Get(obj->GetName());
                if (otherCanvas) {
                    TList* primitives = otherCanvas->GetListOfPrimitives();
                    TIter primIter(primitives);
                    TObject* prim;
                    while ((prim = primIter())) {
                        if (prim->InheritsFrom("TH1")) {
                            TH1* histPrim = (TH1*)prim;
                            TH1* existingHist = (TH1*)canvas->GetPrimitive(histPrim->GetName());
                            if (existingHist) {
                                existingHist->Add(histPrim);
                            } else {
                                canvas->cd();
                                histPrim->Draw("SAME");
                            }
                        }
                    }
                }
            }

            outputFile->cd();
            canvas->Write();
            delete canvas;
        }
    }

    // Chiudi tutti i file
    for (auto* f : files) {
        if (f) f->Close();
    }
    outputFile->Close();

    std::cout << "Histograms and canvases merged successfully into " << outputFileName << std::endl;
}

int combine() {
    // Trova i file automaticamente con un pattern
    std::vector<std::string> inputFiles;
    for (int i = 11; i <= 54; ++i) {
        inputFiles.push_back("out" + std::to_string(i) + ".histAbs.root");
    }
    const char* outputFileName = "A11.combinedhistAbs(dRICHEta).root";

    mergeHisto(inputFiles, outputFileName);

    return 0;
}
