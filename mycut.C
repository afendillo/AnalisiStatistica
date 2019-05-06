#define NSig 10000
#define NBkg 1000000
#define bin 1000
#define rho 0.4
using namespace std;
void Reset(){
	//Recupera la lista di tutti i canvas e li cancella
	TSeqCollection* canvases = gROOT->GetListOfCanvases();
	TIter next(canvases);
	while(TCanvas *c = (TCanvas*)next())
	{
		delete c;
	}
	//Cancella tutti gli oggetti creati
	//Necessario per evitare "warning <TROOT::Append>: Replacing existing TH1: <Name> (Potential memory leak)"
	gROOT->DeleteAll();
	return;
}
string ToString(double n , int precision){
    string i=to_string((int)n);
    string d;
    double decimal=(n-(int)n)*10;
    for(int i=0; i<precision ; i++){
        d+=to_string((int)(decimal));
        decimal=(decimal-(int)decimal)*10;
    }
    return i+"."+d;
}
void mycut(int stamp = 0)
{
	Reset();
	//Informazioni statistiche da stampare
	//gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);
    double PlotRun;

    TFile* MyFile = new TFile("MVA.root");

    if ( MyFile->IsOpen() ) cout<<"File opened successfully\n";
    else {
            cout<<"Error opening MVA.root\n"; 
            return;
        }

    auto eve = (TNtuple *)((MyFile->GetKey("MVA"))->ReadObj());

    TFile* MyCut = new TFile("mycut.root");
    
    if ( MyCut->IsOpen() ) cout<<"Cut opened successfully\n";
    else {
            cout<<"Error opening mycut.root\n"; 
            return;
        }

    auto cutg = (TCutG *)((MyCut->GetKey("CUTG"))->ReadObj());
    cutg->SetVarX("x");
    cutg->SetVarY("y");

    //per poter plottare gli eventi con 4 colori diversi definisco 4 TH2D
    
    TH2D *Sig2D = new TH2D ("Sig2D" , "Taglio Grafico" , bin , -20 , 20 , bin , -20 , 20);
	Sig2D->GetXaxis()->SetTitle("X");Sig2D->GetYaxis()->SetTitle("Y");Sig2D->SetMarkerColor(kGreen);Sig2D->SetMarkerSize(2);Sig2D->SetLineColor(kGreen);

    //eventi background-background -> giallo
    TH2D *Bkg2D = new TH2D ("Bkg2D" , "Bkg2D" , bin , -20 , 20 , bin , -20 , 20);
	Bkg2D->GetXaxis()->SetTitle("X");Bkg2D->GetYaxis()->SetTitle("Y");Bkg2D->SetMarkerColor(kYellow);Bkg2D->SetMarkerSize(2);Bkg2D->SetLineColor(kYellow);

    //eventi background-segnale -> rosso
    TH2D *Pur2D = new TH2D ("Pur2D" , "Pur2D" , bin , -20 , 20 , bin , -20 , 20);
	Pur2D->GetXaxis()->SetTitle("X");Pur2D->GetYaxis()->SetTitle("Y");Pur2D->SetMarkerColor(kRed);Pur2D->SetMarkerSize(2);Pur2D->SetLineColor(kRed);

    //eventi segnale-background -> blu
    TH2D *Eff2D = new TH2D ("Eff2D" , "Eff2D" , bin , -20 , 20 , bin , -20 , 20);
	Eff2D->GetXaxis()->SetTitle("X");Eff2D->GetYaxis()->SetTitle("Y");Eff2D->SetMarkerColor(kBlue);Eff2D->SetMarkerSize(2);Eff2D->SetLineColor(kBlue);

    //disegno sullo stesso canvas gli eventi, divisi in modo opportuno rispetto al taglio
    TCanvas *c2 = new TCanvas();

    eve->Draw("y:x >> Sig2D","s>0 && !CUTG", "");
    eve->Draw("y:x >> Eff2D", "s>0 && CUTG", "SAME");
    eve->Draw("y:x >> Bkg2D", "s<1 && CUTG", "SAME");
    eve->Draw("y:x >> Pur2D", "s<1 && !CUTG", "SAME");
    MyCut->Draw("same");
    
    //calcolo purezza ed efficienza del taglio

    double reiezione, Sig, purezza, efficienza, SS, SB, BS, BB;

   //calcolo purezza ed efficienza del taglio

    SS = Sig2D->GetEntries();
    SB = Eff2D->GetEntries();
    BS = Pur2D->GetEntries();
    BB = Bkg2D->GetEntries();

    purezza = 100*(SS)/(SS+BS);
    efficienza = 100*(SS)/(SS+SB);
    Sig=SS/sqrt(BS+SS);
    reiezione=100*(BB)/(BB+BS);

    TLegend* legend = new TLegend(0.1,0.7,0.45,0.9);
    legend->SetTextSize(0.026); 
    legend->SetTextFont(42);

    legend->SetHeader(("Taglio Grafico #varepsilon: "+ToString(efficienza , 1)+"% p: "+ToString(purezza , 1)+"% r: "+ToString(reiezione , 1)+"%").c_str(),"C");

    legend->AddEntry(Sig2D , "Segnale-Segnale","lp");legend->AddEntry(Eff2D , "Segnale-Fondo","lp");
    legend->AddEntry(Pur2D , "Fondo-Segnale","lp");legend->AddEntry(Bkg2D , "Fondo-Fondo","lp");
    legend->Draw();

    cout<<"##################################################################################\n";
    
    cout<<"Segnale-Segnale: "<<SS<<endl;
    cout<<"Segnale-Fondo: "<<SB<<endl;
    cout<<"Fondo-Fondo: "<<BB<<endl;
    cout<<"Fondo-Segnale: "<<BS<<endl;

    cout <<"Purezza del Segnale Taglio Grafico: "<<purezza<<"%"<<endl;
    cout <<"Efficienza del Segnale Taglio Grafico: "<<efficienza<<"%"<<endl;
    cout <<"Reiezione del fondo Taglio Grafico: "<<reiezione<<"%"<<endl;
    cout <<"Significatività del taglio Taglio Grafico: "<<Sig<<endl;
    cout<<"##################################################################################\n";

    TH1D *EventiX = new TH1D ("EventiX" , "Taglio Grafico Proiezione x" , bin , -10 , 10);
	EventiX->GetXaxis()->SetTitle("X");EventiX->GetYaxis()->SetTitle("Conteggi");

    TH1D *EventiY = new TH1D ("EventiY" , "Taglio Grafico Proiezione y" , bin , -10 , 10);
	EventiY->GetXaxis()->SetTitle("Y");EventiY->GetYaxis()->SetTitle("Conteggi");

    TCanvas *c3 = new TCanvas();
    eve->Draw("x >> EventiX", "!CUTG" , "");

    TCanvas *c4 = new TCanvas();
    eve->Draw("y >> EventiY", "!CUTG" , "");

    if(stamp==1){
        int Cartella1= system("mkdir -p MVA");
        c2->SaveAs("MVA/CutGrafico.png");
        c2->SaveAs("MVA/CutGrafico.root");
        c2->SaveAs("MVA/CutGrafico.pdf");

        c3->SaveAs("MVA/CutGraficoProiezioneX.png");
        c3->SaveAs("MVA/CutGraficoProiezioneX.root");
        c3->SaveAs("MVA/CutGraficoProiezioneX.pdf");

        c4->SaveAs("MVA/CutGraficoProiezioneY.png");
        c4->SaveAs("MVA/CutGraficoProiezioneY.root");
        c4->SaveAs("MVA/CutGraficoProiezioneY.pdf");
    }
 
    return;
}