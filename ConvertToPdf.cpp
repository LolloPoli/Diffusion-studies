#include <TFile.h>
#include <TCanvas.h>
#include <TKey.h>
#include <TH1.h>
#include <TGraph.h>
#include <iostream>
#include <vector>
#include <string>

void ConvertToPDF(const char* filename, const std::vector<int>& indices_to_save) {
    // Apri il file ROOT
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Impossibile aprire il file!" << std::endl;
        return;
    }

    // Ottieni la lista di chiavi nel file
    TList *list_of_keys = file->GetListOfKeys();
    if (!list_of_keys) {
        std::cerr << "Il file non contiene oggetti validi!" << std::endl;
        return;
    }

    // Salva gli oggetti selezionati in un PDF
    TCanvas *c = new TCanvas("c", "Canvas", 800, 600);
    c->Print("output.pdf["); // Inizio del file multipagina

    // Loop sugli indici forniti
    for (int index : indices_to_save) {
        if (index < 0 || index >= list_of_keys->GetEntries()) {
            std::cerr << "Indice fuori intervallo: " << index << std::endl;
            continue;
        }

        // Ottieni l'oggetto dalla lista
        TKey *key = (TKey*)list_of_keys->At(index);
        if (!key) continue;

        TObject *obj = key->ReadObj();
        if (obj->InheritsFrom("TCanvas")) {  // Gestione delle TCanvas
            TCanvas *canvas = (TCanvas*)obj;
            canvas->cd();  // Passa il controllo alla canvas
            canvas->Update();  // Assicurati che il contenuto sia aggiornato
            canvas->Print("output.pdf"); // Salva direttamente la canvas nel PDF
            //std::cout << "Salvato (TCanvas): " << key->GetName() << std::endl;
        } else if (obj->InheritsFrom("TH2")) {  // Gestione degli istogrammi 2D
            obj->Draw("COLZ");
            c->Print("output.pdf");
            //std::cout << "Salvato (TH2): " << key->GetName() << std::endl;
        } else if (obj->InheritsFrom("TH1")) {  // Gestione degli istogrammi 1D
            obj->Draw();
            c->Print("output.pdf");
            //std::cout << "Salvato (TH1): " << key->GetName() << std::endl;
        } else if (obj->InheritsFrom("TGraph")) {  // Gestione dei grafici
            obj->Draw("APL");
            c->Print("output.pdf");
            //std::cout << "Salvato (TGraph): " << key->GetName() << std::endl;
        } else {
            std::cout << "Oggetto ignorato (non grafico): " << key->GetName() << std::endl;
        }

    }

    c->Print("output.pdf]"); // Fine del file multipagina
    file->Close();
    delete c;

    std::cout << "Conversione completata!" << std::endl;
}

int convert() {
    // Nome del file ROOT
    const char* filename = "A10.combinedhistAbs(FullEta).root";

    // Lista di indici degli oggetti da salvare (modifica questa lista)
    std::set<int> exclude_indices = {0, 1, 2, 18, 24, 25, 26, 27, 35, 36, 37, 38, 39, 42, 43, 44, 45, 51, 52, 53, 54, 55};
    //std::set<int> exclude_indices = {111};

    // Lista di indici degli oggetti da salvare
    std::vector<int> indices_to_save;
    for (int i = 0; i <= 64; ++i) {
        if (exclude_indices.count(i)) continue; // Salta se l'indice è nell'insieme
        indices_to_save.push_back(i);
    }

    // Converti in PDF
    ConvertToPDF(filename, indices_to_save);

    return 0;
}
