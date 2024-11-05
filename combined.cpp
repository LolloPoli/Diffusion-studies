#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>

void mergeHisto(const char* inputFile1, const char* inputFile2, const char* inputFile3, const char* inputFile4, 
const char* inputFile5, const char* inputFile6, const char* inputFile7, const char* inputFile8, 
const char* inputFile9, const char* inputFile10, const char* inputFile11, const char* inputFile12,
const char* inputFile13, const char* inputFile14, const char* inputFile15, const char* inputFile16, 
const char* inputFile17, const char* inputFile18, const char* inputFile19, const char* inputFile20, 
const char* inputFile21, const char* inputFile22, const char* inputFile23, const char* inputFile24,
const char* inputFile25, const char* inputFile26, const char* inputFile27,const char* inputFile28, 
const char* inputFile29, const char* inputFile30, const char* outputFileName) {
    // Open the first input file
    TFile *file1 = TFile::Open(inputFile1);
    if (!file1 || file1->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile1 << std::endl;
        return;
    }

    // Open the second input file
    TFile *file2 = TFile::Open(inputFile2);
    if (!file2 || file2->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile2 << std::endl;
        file1->Close();
        return;
    }

    // Open the third input file
    TFile *file3 = TFile::Open(inputFile3);
    if (!file3 || file3->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile3 << std::endl;
        file1->Close();
        file2->Close();
        return;
    }

    // Open the fourth input file
    TFile *file4 = TFile::Open(inputFile4);
    if (!file4 || file4->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile4 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        return;
    }

    // Open the fifth input file
    TFile *file5 = TFile::Open(inputFile5);
    if (!file5 || file5->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile5 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        return;
    }

    // Open the fifth input file
    TFile *file6 = TFile::Open(inputFile6);
    if (!file6 || file6->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile6 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        return;
    }

     // Open the fifth input file
    TFile *file7 = TFile::Open(inputFile7);
    if (!file7 || file7->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile7 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        return;
    }

    TFile *file8 = TFile::Open(inputFile8);
    if (!file8 || file8->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile8 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        return;
    }

    TFile *file9 = TFile::Open(inputFile9);
    if (!file9 || file9->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile9 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        return;
    }

    TFile *file10 = TFile::Open(inputFile10);
    if (!file10 || file10->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile10 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        return;
    }

    TFile *file11 = TFile::Open(inputFile11);
    if (!file11 || file11->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile11 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        return;
    }

    TFile *file12 = TFile::Open(inputFile12);
    if (!file12 || file12->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile12 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        return;
    }

    TFile *file13 = TFile::Open(inputFile13);
    if (!file13 || file13->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile13 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        return;
    }

    TFile *file14 = TFile::Open(inputFile14);
    if (!file14 || file14->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile14 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        return;
    }
    
    TFile *file15 = TFile::Open(inputFile15);
    if (!file15 || file15->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile15 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        return;
    }

    TFile *file16 = TFile::Open(inputFile16);
    if (!file16 || file16->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile16 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        file15->Close();
        return;
    }

    TFile *file17 = TFile::Open(inputFile17);
    if (!file17 || file17->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile17 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        file15->Close();
        file16->Close();
        return;
    }

    TFile *file18 = TFile::Open(inputFile18);
    if (!file18 || file18->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile18 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        file15->Close();
        file16->Close();
        file17->Close();
        return;
    }

    TFile *file19 = TFile::Open(inputFile19);
    if (!file19 || file19->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile19 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        file15->Close();
        file16->Close();
        file17->Close();
        file18->Close();
        return;
    }

    TFile *file20 = TFile::Open(inputFile20);
    if (!file20 || file20->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile20 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        file15->Close();
        file16->Close();
        file17->Close();
        file18->Close();
        file19->Close();
        return;
    }

    TFile *file21 = TFile::Open(inputFile21);
    if (!file21 || file21->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile21 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        file15->Close();
        file16->Close();
        file17->Close();
        file18->Close();
        file19->Close();
        file20->Close();
        return;
    }

    TFile *file22 = TFile::Open(inputFile22);
    if (!file22 || file22->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile22 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        file15->Close();
        file16->Close();
        file17->Close();
        file18->Close();
        file19->Close();
        file20->Close();
        file21->Close();
        return;
    }

    TFile *file23 = TFile::Open(inputFile23);
    if (!file23 || file23->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile23 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        file15->Close();
        file16->Close();
        file17->Close();
        file18->Close();
        file19->Close();
        file20->Close();
        file21->Close();
        file22->Close();
        return;
    }

    TFile *file24 = TFile::Open(inputFile24);
    if (!file24 || file24->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile24 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        file15->Close();
        file16->Close();
        file17->Close();
        file18->Close();
        file19->Close();
        file20->Close();
        file21->Close();
        file22->Close();
        file23->Close();
        return;
    }

    TFile *file25 = TFile::Open(inputFile25);
    if (!file25 || file25->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile25 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        file15->Close();
        file16->Close();
        file17->Close();
        file18->Close();
        file19->Close();
        file20->Close();
        file21->Close();
        file22->Close();
        file23->Close();
        file24->Close();
        return;
    }

    TFile *file26 = TFile::Open(inputFile26);
    if (!file26 || file26->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile26 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        file15->Close();
        file16->Close();
        file17->Close();
        file18->Close();
        file19->Close();
        file20->Close();
        file21->Close();
        file22->Close();
        file23->Close();
        file24->Close();
        file25->Close();
        return;
    }

    TFile *file27 = TFile::Open(inputFile27);
    if (!file27 || file27->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile27 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        file15->Close();
        file16->Close();
        file17->Close();
        file18->Close();
        file19->Close();
        file20->Close();
        file21->Close();
        file22->Close();
        file23->Close();
        file24->Close();
        file25->Close();
        file26->Close();
        return;
    }

    TFile *file28 = TFile::Open(inputFile28);
    if (!file28 || file28->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile28 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        file15->Close();
        file16->Close();
        file17->Close();
        file18->Close();
        file19->Close();
        file20->Close();
        file21->Close();
        file22->Close();
        file23->Close();
        file24->Close();
        file25->Close();
        file26->Close();
        file27->Close();
        return;
    }

    TFile *file29 = TFile::Open(inputFile29);
    if (!file29 || file29->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile29 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        file15->Close();
        file16->Close();
        file17->Close();
        file18->Close();
        file19->Close();
        file20->Close();
        file21->Close();
        file22->Close();
        file23->Close();
        file24->Close();
        file25->Close();
        file26->Close();
        file27->Close();
        file28->Close();
        return;
    }

    TFile *file30 = TFile::Open(inputFile30);
    if (!file30 || file30->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputFile30 << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        file15->Close();
        file16->Close();
        file17->Close();
        file18->Close();
        file19->Close();
        file20->Close();
        file21->Close();
        file22->Close();
        file23->Close();
        file24->Close();
        file25->Close();
        file26->Close();
        file27->Close();
        file28->Close();
        file29->Close();
        return;
    }

    // Open the output file in update mode
    TFile *outputFile = TFile::Open(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Cannot open or create output file " << outputFileName << std::endl;
        file1->Close();
        file2->Close();
        file3->Close();
        file4->Close();
        file5->Close();
        file6->Close();
        file7->Close();
        file8->Close();
        file9->Close();
        file10->Close();
        file11->Close();
        file12->Close();
        file13->Close();
        file14->Close();
        file15->Close();
        file16->Close();
        file17->Close();
        file18->Close();
        file19->Close();
        file20->Close();
        file21->Close();
        file22->Close();
        file23->Close();
        file24->Close();
        file25->Close();
        file26->Close();
        file27->Close();
        file28->Close();
        file29->Close();
        file30->Close();
        return;
    }   
    

    // Merge canvases
    TList *canvasList = file1->GetListOfKeys();
    TIter nextCanvas(canvasList);
    TKey *keyCanvas;
    while ((keyCanvas = (TKey*)nextCanvas())) {
        TObject *objCanvas = keyCanvas->ReadObj();
        if (objCanvas->InheritsFrom("TCanvas")) {
            TCanvas *canvas1 = (TCanvas*)objCanvas;
            TCanvas *canvas2 = (TCanvas*)file2->Get(canvas1->GetName());
            TCanvas *canvas3 = (TCanvas*)file3->Get(canvas1->GetName());
            TCanvas *canvas4 = (TCanvas*)file4->Get(canvas1->GetName());
            TCanvas *canvas5 = (TCanvas*)file5->Get(canvas1->GetName());
            TCanvas *canvas6 = (TCanvas*)file6->Get(canvas1->GetName());
            TCanvas *canvas7 = (TCanvas*)file7->Get(canvas1->GetName());
            TCanvas *canvas8 = (TCanvas*)file8->Get(canvas1->GetName());
            TCanvas *canvas9 = (TCanvas*)file9->Get(canvas1->GetName());
            TCanvas *canvas10 = (TCanvas*)file10->Get(canvas1->GetName());
            TCanvas *canvas11 = (TCanvas*)file11->Get(canvas1->GetName());
            TCanvas *canvas12 = (TCanvas*)file12->Get(canvas1->GetName());
            TCanvas *canvas13 = (TCanvas*)file13->Get(canvas1->GetName());
            TCanvas *canvas14 = (TCanvas*)file14->Get(canvas1->GetName());
            TCanvas *canvas15 = (TCanvas*)file15->Get(canvas1->GetName());
            TCanvas *canvas16 = (TCanvas*)file16->Get(canvas1->GetName());
            TCanvas *canvas17 = (TCanvas*)file17->Get(canvas1->GetName());
            TCanvas *canvas18 = (TCanvas*)file18->Get(canvas1->GetName());
            TCanvas *canvas19 = (TCanvas*)file19->Get(canvas1->GetName());
            TCanvas *canvas20 = (TCanvas*)file20->Get(canvas1->GetName());
            TCanvas *canvas21 = (TCanvas*)file21->Get(canvas1->GetName());
            TCanvas *canvas22 = (TCanvas*)file22->Get(canvas1->GetName());
            TCanvas *canvas23 = (TCanvas*)file23->Get(canvas1->GetName());
            TCanvas *canvas24 = (TCanvas*)file24->Get(canvas1->GetName());
            TCanvas *canvas25 = (TCanvas*)file25->Get(canvas1->GetName());
            TCanvas *canvas26 = (TCanvas*)file26->Get(canvas1->GetName());
            TCanvas *canvas27 = (TCanvas*)file27->Get(canvas1->GetName());
            TCanvas *canvas28 = (TCanvas*)file28->Get(canvas1->GetName());
            TCanvas *canvas29 = (TCanvas*)file29->Get(canvas1->GetName());
            TCanvas *canvas30 = (TCanvas*)file30->Get(canvas1->GetName());
            
            if (canvas2) {
                TList *listPrimitives2 = canvas2->GetListOfPrimitives();
                TIter nextPrim2(listPrimitives2);
                TObject *objPrim2;
                while ((objPrim2 = nextPrim2())) {
                    if (objPrim2->InheritsFrom("TH1")) {
                        TH1 *hist2 = (TH1*)objPrim2;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist2->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist2);
                        } else {
                            canvas1->cd();
                            hist2->Draw("SAME");
                        }
                    }
                }
            }
            
            if (canvas3) {
                TList *listPrimitives3 = canvas3->GetListOfPrimitives();
                TIter nextPrim3(listPrimitives3);
                TObject *objPrim3;
                while ((objPrim3 = nextPrim3())) {
                    if (objPrim3->InheritsFrom("TH1")) {
                        TH1 *hist3 = (TH1*)objPrim3;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist3->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist3);
                        } else {
                            canvas1->cd();
                            hist3->Draw("SAME");
                        }
                    }
                }
            }
            
            if (canvas4) {
                TList *listPrimitives4 = canvas4->GetListOfPrimitives();
                TIter nextPrim4(listPrimitives4);
                TObject *objPrim4;
                while ((objPrim4 = nextPrim4())) {
                    if (objPrim4->InheritsFrom("TH1")) {
                        TH1 *hist4 = (TH1*)objPrim4;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist4->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist4);
                        } else {
                            canvas1->cd();
                            hist4->Draw("SAME");
                        }
                    }
                }
            }
            
            if (canvas5) {
                TList *listPrimitives5 = canvas5->GetListOfPrimitives();
                TIter nextPrim5(listPrimitives5);
                TObject *objPrim5;
                while ((objPrim5 = nextPrim5())) {
                    if (objPrim5->InheritsFrom("TH1")) {
                        TH1 *hist5 = (TH1*)objPrim5;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist5->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist5);
                        } else {
                            canvas1->cd();
                            hist5->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas6) {
                TList *listPrimitives6 = canvas6->GetListOfPrimitives();
                TIter nextPrim6(listPrimitives6);
                TObject *objPrim6;
                while ((objPrim6 = nextPrim6())) {
                    if (objPrim6->InheritsFrom("TH1")) {
                        TH1 *hist6 = (TH1*)objPrim6;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist6->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist6);
                        } else {
                            canvas1->cd();
                            hist6->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas7) {
                TList *listPrimitives7 = canvas7->GetListOfPrimitives();
                TIter nextPrim7(listPrimitives7);
                TObject *objPrim7;
                while ((objPrim7 = nextPrim7())) {
                    if (objPrim7->InheritsFrom("TH1")) {
                        TH1 *hist7 = (TH1*)objPrim7;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist7->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist7);
                        } else {
                            canvas1->cd();
                            hist7->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas8) {
                TList *listPrimitives8 = canvas8->GetListOfPrimitives();
                TIter nextPrim8(listPrimitives8);
                TObject *objPrim8;
                while ((objPrim8 = nextPrim8())) {
                    if (objPrim8->InheritsFrom("TH1")) {
                        TH1 *hist8 = (TH1*)objPrim8;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist8->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist8);
                        } else {
                            canvas1->cd();
                            hist8->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas9) {
                TList *listPrimitives9 = canvas9->GetListOfPrimitives();
                TIter nextPrim9(listPrimitives9);
                TObject *objPrim9;
                while ((objPrim9 = nextPrim9())) {
                    if (objPrim9->InheritsFrom("TH1")) {
                        TH1 *hist9 = (TH1*)objPrim9;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist9->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist9);
                        } else {
                            canvas1->cd();
                            hist9->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas10) {
                TList *listPrimitives10 = canvas10->GetListOfPrimitives();
                TIter nextPrim10(listPrimitives10);
                TObject *objPrim10;
                while ((objPrim10 = nextPrim10())) {
                    if (objPrim10->InheritsFrom("TH1")) {
                        TH1 *hist10 = (TH1*)objPrim10;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist10->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist10);
                        } else {
                            canvas1->cd();
                            hist10->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas11) {
                TList *listPrimitives11 = canvas11->GetListOfPrimitives();
                TIter nextPrim11(listPrimitives11);
                TObject *objPrim11;
                while ((objPrim11 = nextPrim11())) {
                    if (objPrim11->InheritsFrom("TH1")) {
                        TH1 *hist11 = (TH1*)objPrim11;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist11->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist11);
                        } else {
                            canvas1->cd();
                            hist11->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas12) {
                TList *listPrimitives12 = canvas12->GetListOfPrimitives();
                TIter nextPrim12(listPrimitives12);
                TObject *objPrim12;
                while ((objPrim12 = nextPrim12())) {
                    if (objPrim12->InheritsFrom("TH1")) {
                        TH1 *hist12 = (TH1*)objPrim12;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist12->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist12);
                        } else {
                            canvas1->cd();
                            hist12->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas13) {
                TList *listPrimitives13 = canvas13->GetListOfPrimitives();
                TIter nextPrim13(listPrimitives13);
                TObject *objPrim13;
                while ((objPrim13 = nextPrim13())) {
                    if (objPrim13->InheritsFrom("TH1")) {
                        TH1 *hist13 = (TH1*)objPrim13;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist13->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist13);
                        } else {
                            canvas1->cd();
                            hist13->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas14) {
                TList *listPrimitives14 = canvas14->GetListOfPrimitives();
                TIter nextPrim14(listPrimitives14);
                TObject *objPrim14;
                while ((objPrim14 = nextPrim14())) {
                    if (objPrim14->InheritsFrom("TH1")) {
                        TH1 *hist14 = (TH1*)objPrim14;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist14->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist14);
                        } else {
                            canvas1->cd();
                            hist14->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas15) {
                TList *listPrimitives15 = canvas15->GetListOfPrimitives();
                TIter nextPrim15(listPrimitives15);
                TObject *objPrim15;
                while ((objPrim15 = nextPrim15())) {
                    if (objPrim15->InheritsFrom("TH1")) {
                        TH1 *hist15 = (TH1*)objPrim15;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist15->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist15);
                        } else {
                            canvas1->cd();
                            hist15->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas16) {
                TList *listPrimitives16 = canvas16->GetListOfPrimitives();
                TIter nextPrim16(listPrimitives16);
                TObject *objPrim16;
                while ((objPrim16 = nextPrim16())) {
                    if (objPrim16->InheritsFrom("TH1")) {
                        TH1 *hist16 = (TH1*)objPrim16;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist16->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist16);
                        } else {
                            canvas1->cd();
                            hist16->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas17) {
                TList *listPrimitives17 = canvas17->GetListOfPrimitives();
                TIter nextPrim17(listPrimitives17);
                TObject *objPrim17;
                while ((objPrim17 = nextPrim17())) {
                    if (objPrim17->InheritsFrom("TH1")) {
                        TH1 *hist17 = (TH1*)objPrim17;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist17->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist17);
                        } else {
                            canvas1->cd();
                            hist17->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas18) {
                TList *listPrimitives18 = canvas18->GetListOfPrimitives();
                TIter nextPrim18(listPrimitives18);
                TObject *objPrim18;
                while ((objPrim18 = nextPrim18())) {
                    if (objPrim18->InheritsFrom("TH1")) {
                        TH1 *hist18 = (TH1*)objPrim18;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist18->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist18);
                        } else {
                            canvas1->cd();
                            hist18->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas19) {
                TList *listPrimitives19 = canvas19->GetListOfPrimitives();
                TIter nextPrim19(listPrimitives19);
                TObject *objPrim19;
                while ((objPrim19 = nextPrim19())) {
                    if (objPrim19->InheritsFrom("TH1")) {
                        TH1 *hist19 = (TH1*)objPrim19;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist19->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist19);
                        } else {
                            canvas1->cd();
                            hist19->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas20) {
                TList *listPrimitives20 = canvas20->GetListOfPrimitives();
                TIter nextPrim20(listPrimitives20);
                TObject *objPrim20;
                while ((objPrim20 = nextPrim20())) {
                    if (objPrim20->InheritsFrom("TH1")) {
                        TH1 *hist20 = (TH1*)objPrim20;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist20->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist20);
                        } else {
                            canvas1->cd();
                            hist20->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas21) {
                TList *listPrimitives21 = canvas21->GetListOfPrimitives();
                TIter nextPrim21(listPrimitives21);
                TObject *objPrim21;
                while ((objPrim21 = nextPrim21())) {
                    if (objPrim21->InheritsFrom("TH1")) {
                        TH1 *hist21 = (TH1*)objPrim21;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist21->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist21);
                        } else {
                            canvas1->cd();
                            hist21->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas22) {
                TList *listPrimitives22 = canvas22->GetListOfPrimitives();
                TIter nextPrim22(listPrimitives22);
                TObject *objPrim22;
                while ((objPrim22 = nextPrim22())) {
                    if (objPrim22->InheritsFrom("TH1")) {
                        TH1 *hist22 = (TH1*)objPrim22;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist22->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist22);
                        } else {
                            canvas1->cd();
                            hist22->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas23) {
                TList *listPrimitives23 = canvas23->GetListOfPrimitives();
                TIter nextPrim23(listPrimitives23);
                TObject *objPrim23;
                while ((objPrim23 = nextPrim23())) {
                    if (objPrim23->InheritsFrom("TH1")) {
                        TH1 *hist23 = (TH1*)objPrim23;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist23->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist23);
                        } else {
                            canvas1->cd();
                            hist23->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas24) {
                TList *listPrimitives24 = canvas24->GetListOfPrimitives();
                TIter nextPrim24(listPrimitives24);
                TObject *objPrim24;
                while ((objPrim24 = nextPrim24())) {
                    if (objPrim24->InheritsFrom("TH1")) {
                        TH1 *hist24 = (TH1*)objPrim24;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist24->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist24);
                        } else {
                            canvas1->cd();
                            hist24->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas25) {
                TList *listPrimitives25 = canvas25->GetListOfPrimitives();
                TIter nextPrim25(listPrimitives25);
                TObject *objPrim25;
                while ((objPrim25 = nextPrim25())) {
                    if (objPrim25->InheritsFrom("TH1")) {
                        TH1 *hist25 = (TH1*)objPrim25;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist25->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist25);
                        } else {
                            canvas1->cd();
                            hist25->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas26) {
                TList *listPrimitives26 = canvas26->GetListOfPrimitives();
                TIter nextPrim26(listPrimitives26);
                TObject *objPrim26;
                while ((objPrim26 = nextPrim26())) {
                    if (objPrim26->InheritsFrom("TH1")) {
                        TH1 *hist26 = (TH1*)objPrim26;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist26->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist26);
                        } else {
                            canvas1->cd();
                            hist26->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas27) {
                TList *listPrimitives27 = canvas27->GetListOfPrimitives();
                TIter nextPrim27(listPrimitives27);
                TObject *objPrim27;
                while ((objPrim27 = nextPrim27())) {
                    if (objPrim27->InheritsFrom("TH1")) {
                        TH1 *hist27 = (TH1*)objPrim27;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist27->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist27);
                        } else {
                            canvas1->cd();
                            hist27->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas28) {
                TList *listPrimitives28 = canvas28->GetListOfPrimitives();
                TIter nextPrim28(listPrimitives28);
                TObject *objPrim28;
                while ((objPrim28 = nextPrim28())) {
                    if (objPrim28->InheritsFrom("TH1")) {
                        TH1 *hist28 = (TH1*)objPrim28;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist28->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist28);
                        } else {
                            canvas1->cd();
                            hist28->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas29) {
                TList *listPrimitives29 = canvas29->GetListOfPrimitives();
                TIter nextPrim29(listPrimitives29);
                TObject *objPrim29;
                while ((objPrim29 = nextPrim29())) {
                    if (objPrim29->InheritsFrom("TH1")) {
                        TH1 *hist29 = (TH1*)objPrim29;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist29->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist29);
                        } else {
                            canvas1->cd();
                            hist29->Draw("SAME");
                        }
                    }
                }
            }

            if (canvas30) {
                TList *listPrimitives30 = canvas30->GetListOfPrimitives();
                TIter nextPrim30(listPrimitives30);
                TObject *objPrim30;
                while ((objPrim30 = nextPrim30())) {
                    if (objPrim30->InheritsFrom("TH1")) {
                        TH1 *hist30 = (TH1*)objPrim30;
                        TH1 *existingHist1 = (TH1*)canvas1->GetPrimitive(hist30->GetName());
                        if (existingHist1) {
                            existingHist1->Add(hist30);
                        } else {
                            canvas1->cd();
                            hist30->Draw("SAME");
                        }
                    }
                }
            }
            
            outputFile->cd();
            canvas1->Write(canvas1->GetName(), TObject::kOverwrite);
        }
    }


    // Merge histograms
    TList *keyList = file1->GetListOfKeys();
    TIter next(keyList);
    TKey *key;
    while ((key = (TKey*)next())) {
        TObject *obj = key->ReadObj();
        if (obj->InheritsFrom("TH1")) {
            TH1 *hist1 = (TH1*)obj;
            TH1 *hist2 = (TH1*)file2->Get(hist1->GetName());
            TH1 *hist3 = (TH1*)file3->Get(hist1->GetName());
            TH1 *hist4 = (TH1*)file4->Get(hist1->GetName());
            TH1 *hist5 = (TH1*)file5->Get(hist1->GetName());
            TH1 *hist6 = (TH1*)file6->Get(hist1->GetName());
            TH1 *hist7 = (TH1*)file7->Get(hist1->GetName());
            TH1 *hist8 = (TH1*)file8->Get(hist1->GetName());
            TH1 *hist9 = (TH1*)file9->Get(hist1->GetName());
            TH1 *hist10 = (TH1*)file10->Get(hist1->GetName());
            TH1 *hist11 = (TH1*)file11->Get(hist1->GetName());
            TH1 *hist12 = (TH1*)file12->Get(hist1->GetName());
            TH1 *hist13 = (TH1*)file13->Get(hist1->GetName());
            TH1 *hist14 = (TH1*)file14->Get(hist1->GetName());
            TH1 *hist15 = (TH1*)file15->Get(hist1->GetName());
            TH1 *hist16 = (TH1*)file16->Get(hist1->GetName());
            TH1 *hist17 = (TH1*)file17->Get(hist1->GetName());
            TH1 *hist18 = (TH1*)file18->Get(hist1->GetName());
            TH1 *hist19 = (TH1*)file19->Get(hist1->GetName());
            TH1 *hist20 = (TH1*)file20->Get(hist1->GetName());
            TH1 *hist21 = (TH1*)file21->Get(hist1->GetName());
            TH1 *hist22 = (TH1*)file22->Get(hist1->GetName());
            TH1 *hist23 = (TH1*)file23->Get(hist1->GetName());
            TH1 *hist24 = (TH1*)file24->Get(hist1->GetName());
            TH1 *hist25 = (TH1*)file25->Get(hist1->GetName());
            TH1 *hist26 = (TH1*)file26->Get(hist1->GetName());
            TH1 *hist27 = (TH1*)file27->Get(hist1->GetName());
            TH1 *hist28 = (TH1*)file28->Get(hist1->GetName());
            TH1 *hist29 = (TH1*)file29->Get(hist1->GetName());
            TH1 *hist30 = (TH1*)file30->Get(hist1->GetName());
            
            if (hist2) {
                hist1->Add(hist2);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile2 << std::endl;
            }
            
            if (hist3) {
                hist1->Add(hist3);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile3 << std::endl;
            }
            
            if (hist4) {
                hist1->Add(hist4);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile4 << std::endl;
            }
            
            if (hist5) {
                hist1->Add(hist5);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile5 << std::endl;
            }

            if (hist6) {
                hist1->Add(hist6);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile6 << std::endl;
            }

            if (hist7) {
                hist1->Add(hist7);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile7 << std::endl;
            }

            if (hist8) {
                hist1->Add(hist8);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile8 << std::endl;
            }

            if (hist9) {
                hist1->Add(hist9);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile9 << std::endl;
            }

            if (hist10) {
                hist1->Add(hist10);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile10 << std::endl;
            }

            if (hist11) {
                hist1->Add(hist11);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile11 << std::endl;
            }

            if (hist12) {
                hist1->Add(hist12);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile12 << std::endl;
            }

            if (hist13) {
                hist1->Add(hist13);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile13 << std::endl;
            }   

            if (hist14) {
                hist1->Add(hist14);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile14 << std::endl;
            } 

            if (hist15) {
                hist1->Add(hist15);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile15 << std::endl;
            }

            if (hist16) {
                hist1->Add(hist16);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile16 << std::endl;
            }  

            if (hist17) {
                hist1->Add(hist17);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile17 << std::endl;
            }

            if (hist18) {
                hist1->Add(hist18);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile18 << std::endl;
            }

            if (hist19) {
                hist1->Add(hist19);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile19 << std::endl;
            }

            if (hist20) {
                hist1->Add(hist20);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile20 << std::endl;
            }  

            if (hist21) {
                hist1->Add(hist21);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile21 << std::endl;
            }

            if (hist22) {
                hist1->Add(hist22);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile22 << std::endl;
            }

            if (hist23) {
                hist1->Add(hist23);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile23 << std::endl;
            } 

            if (hist24) {
                hist1->Add(hist24);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile24 << std::endl;
            } 

            if (hist25) {
                hist1->Add(hist25);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile25 << std::endl;
            } 

            if (hist26) {
                hist1->Add(hist26);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile26 << std::endl;
            } 

            if (hist27) {
                hist1->Add(hist27);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile27 << std::endl;
            } 

            if (hist28) {
                hist1->Add(hist28);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile28 << std::endl;
            } 

            if (hist29) {
                hist1->Add(hist29);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile29 << std::endl;
            } 

            if (hist30) {
                hist1->Add(hist30);
            } else {
                std::cerr << "Warning: Histogram " << hist1->GetName() << " not found in " << inputFile30 << std::endl;
            } 
            
            outputFile->cd();
            hist1->Write(hist1->GetName(), TObject::kOverwrite);
        }
    }

    // Close all files
    file1->Close();
    file2->Close();
    file3->Close();
    file4->Close();
    file5->Close();
    file6->Close();
    file7->Close();
    file8->Close();
    file9->Close();
    file10->Close();
    file11->Close();
    file12->Close();
    file13->Close();
    file14->Close();
    file15->Close();
    file16->Close();
    file17->Close();
    file18->Close();
    file19->Close();
    file20->Close();
    file21->Close();
    file22->Close();
    file23->Close();
    file24->Close();
    file25->Close();
    file26->Close();
    file27->Close();
    file28->Close();
    file29->Close();
    file30->Close();
    outputFile->Close();

    std::cout << "Histograms and canvases merged successfully into " << outputFileName << std::endl;
}

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
/*
int main() {
    const char* inputFile1 = "out1.hist.root";
    const char* inputFile2 = "out2.hist.root";
    const char* inputFile3 = "out3.hist.root";
    const char* inputFile4 = "out4.hist.root";
    const char* inputFile5 = "out5.hist.root";
    const char* outputFileName = "combined.hist.root";

    mergeHisto(inputFile1, inputFile2, inputFile3, inputFile4, inputFile5, outputFileName);
    // mergeHisto("out1.hist.root", "out2.hist.root", "out3.hist.root", "out4.hist.root", "out5.hist.root", "combined.hist.root")
    // mergeHisto("out6.hist.root", "out7.hist.root", "out8.hist.root", "out9.hist.root", "out10.hist.root", "combined.hist2.root")

    return 0;
}

int main2() {
    const char* inputFile1 = "out1.histN.root";
    const char* inputFile2 = "out2.histN.root";
    const char* inputFile3 = "out3.histN.root";
    const char* inputFile4 = "out4.histN.root";
    const char* inputFile5 = "out5.histN.root";
    const char* outputFileName = "combined.histAbs.root";

    mergeHisto(inputFile1, inputFile2, inputFile3, inputFile4, inputFile5, outputFileName);

    return 0;
}
*/
int mainPos() {
    const char* inputFile1 = "out11.hist.root";
    const char* inputFile2 = "out12.hist.root";
    const char* inputFile3 = "out13.hist.root";
    const char* inputFile4 = "out14.hist.root";
    const char* inputFile5 = "out15.hist.root";
    const char* inputFile6 = "out16.hist.root";
    const char* inputFile7 = "out17.hist.root";
    const char* inputFile8 = "out18.hist.root";
    const char* inputFile9 = "out19.hist.root";
    const char* inputFile10 = "out20.hist.root";
    const char* inputFile11 = "out21.hist.root";
    const char* inputFile12 = "out22.hist.root";
    const char* inputFile13 = "out23.hist.root";
    const char* inputFile14 = "out24.hist.root";
    const char* inputFile15 = "out25.hist.root";
    const char* inputFile16 = "out26.hist.root";
    const char* inputFile17 = "out27.hist.root";
    const char* inputFile18 = "out28.hist.root";
    const char* inputFile19 = "out29.hist.root";
    const char* inputFile20 = "out30.hist.root";
    const char* inputFile21 = "out31.hist.root";
    const char* inputFile22 = "out32.hist.root";
    const char* inputFile23 = "out33.hist.root";
    const char* inputFile24 = "out34.hist.root";
    const char* inputFile25 = "out35.hist.root";
    const char* inputFile26 = "out36.hist.root";
    const char* inputFile27 = "out37.hist.root";
    const char* inputFile28 = "out38.hist.root";
    const char* inputFile29 = "out39.hist.root";
    const char* inputFile30 = "out40.hist.root";
    const char* outputFileName = "combined.hist2.root";

    mergeHisto(inputFile1, inputFile2, inputFile3, inputFile4, inputFile5, inputFile6, inputFile7, inputFile8, 
    inputFile9, inputFile10, inputFile11, inputFile12, inputFile13,  inputFile14, inputFile15, inputFile16,
    inputFile17, inputFile18, inputFile19, inputFile20, inputFile21, inputFile22, inputFile23, inputFile24, 
    inputFile25, inputFile26, inputFile27, inputFile28, inputFile29, inputFile30, outputFileName);

    return 0;
}

int mainAbs() {
    const char* inputFile1 = "out11.histAbs.root";
    const char* inputFile2 = "out12.histAbs.root";
    const char* inputFile3 = "out13.histAbs.root";
    const char* inputFile4 = "out14.histAbs.root";
    const char* inputFile5 = "out15.histAbs.root";
    const char* inputFile6 = "out16.histAbs.root";
    const char* inputFile7 = "out17.histAbs.root";
    const char* inputFile8 = "out18.histAbs.root";
    const char* inputFile9 = "out19.histAbs.root";
    const char* inputFile10 = "out20.histAbs.root";
    const char* inputFile11 = "out21.histAbs.root";
    const char* inputFile12 = "out22.histAbs.root";
    const char* inputFile13 = "out23.histAbs.root";
    const char* inputFile14 = "out24.histAbs.root";
    const char* inputFile15 = "out25.histAbs.root";
    const char* inputFile16 = "out26.histAbs.root";
    const char* inputFile17 = "out27.histAbs.root";
    const char* inputFile18 = "out28.histAbs.root";
    const char* inputFile19 = "out29.histAbs.root";
    const char* inputFile20 = "out30.histAbs.root";
    const char* inputFile21 = "out31.histAbs.root";
    const char* inputFile22 = "out32.histAbs.root";
    const char* inputFile23 = "out33.histAbs.root";
    const char* inputFile24 = "out34.histAbs.root";
    const char* inputFile25 = "out35.histAbs.root";
    const char* inputFile26 = "out36.histAbs.root";
    const char* inputFile27 = "out37.histAbs.root";
    const char* inputFile28 = "out38.histAbs.root";
    const char* inputFile29 = "out39.histAbs.root";
    const char* inputFile30 = "out40.histAbs.root";
    const char* outputFileName = "A.combinedhistAbs.root";

    mergeHisto(inputFile1, inputFile2, inputFile3, inputFile4, inputFile5, inputFile6, inputFile7, inputFile8, 
    inputFile9, inputFile10, inputFile11, inputFile12, inputFile13,  inputFile14, inputFile15, inputFile16,
    inputFile17, inputFile18, inputFile19, inputFile20, inputFile21, inputFile22, inputFile23, inputFile24, 
    inputFile25, inputFile26, inputFile27, inputFile28, inputFile29, inputFile30, outputFileName);

    return 0;
}
