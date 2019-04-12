#define NSig 10000
#define NBkg 1000000
#define bin 1000
#define rho 0.4
using namespace std;


void mycut(char stamp='Q')
{
    //Recupera la lista di tutti i canvas e li cancella
	int count=0;
	TSeqCollection* canvases = gROOT->GetListOfCanvases();
	TIter next(canvases);
	while(TCanvas *c = (TCanvas*)next())
	{
		delete c;
		count++;
	}
	if(stamp!='Q') cout<<count<<" Canvases Deleted"<<endl;
	//Cancella tutti gli oggetti creati
	//Necessario per evitare "warning <TROOT::Append>: Replacing existing TH1: <Name> (Potential memory leak)"
	gROOT->DeleteAll();
	
	//Informazioni statistiche da stampare
	//gStyle->SetOptFit(1111);

    TNtuple *eve = new TNtuple("events", "events", "x:y:s"); //ebbene si, usiamo le NTuple

    TRandom3 *Rand = new TRandom3(time(0)); 

    double x, y; int s;

    for (int i = 0; i<NSig; i++) //Genera il segnale
    {
        x = Rand->BreitWigner(0,1);
        y = Rand->BreitWigner(1,1);
        if(abs(x)<20 && abs(y)<20) //per ottenere un grafico leggibile dobbiamo restringere attorno a 0
            eve->Fill(x,y,1);
        else i--;
    }

    for (int i = 0; i<NBkg; i++) //Genera il background
    {
        //per generare una coppia di variabili correlate con correlazione rho 
        //ci serve prima una coppia di variabili scorrelate con varianza unitaria
        x = Rand->Gaus(0,1);
        y = Rand->Gaus(0,1);
        if(abs(x)<20 && abs(y)<20) //solito controllo
        {
            //quindi la combinazione di x
            y = y*sqrt(1-(rho*rho)) + x*rho;
            //sarà correlata con indice rho e manterrà varianza unitaria
        
            //centriamo ora la gaussiana in 4,4
            x+=4;
            y+=4;
            eve->Fill(x,y,0);
        }
        else i--;
    }

    //per poter plottare gli eventi con 4 colori diversi definisco 4 TH2D
    
    //eventi segnale-segnale -> verde
    TH2D *Sig2D = new TH2D ("Sig2D" , "Plot Run" , bin , -20 , 20 , bin , -20 , 20);
	Sig2D->GetXaxis()->SetTitle("X");
	Sig2D->GetYaxis()->SetTitle("Y");
	Sig2D->SetMarkerColor(kGreen);

    //eventi background-background -> giallo
    TH2D *Bkg2D = new TH2D ("Bkg2D" , "Bkg2D" , bin , -20 , 20 , bin , -20 , 20);
	Bkg2D->GetXaxis()->SetTitle("X");
	Bkg2D->GetYaxis()->SetTitle("Y");
	Bkg2D->SetMarkerColor(kYellow);

    //eventi background-segnale -> rosso
    TH2D *Pur2D = new TH2D ("Pur2D" , "Pur2D" , bin , -20 , 20 , bin , -20 , 20);
	Pur2D->GetXaxis()->SetTitle("X");
	Pur2D->GetYaxis()->SetTitle("Y");
	Pur2D->SetMarkerColor(kRed);

    //eventi segnale-background -> blu
    TH2D *Eff2D = new TH2D ("Eff2D" , "Eff2D" , bin , -20 , 20 , bin , -20 , 20);
	Eff2D->GetXaxis()->SetTitle("X");
	Eff2D->GetYaxis()->SetTitle("Y");
	Eff2D->SetMarkerColor(kBlue);

     
   TCutG *cutg = new TCutG("CUTG",14);
   cutg->SetVarX("x");
   cutg->SetVarY("y");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,0.286533,3.88655);
   cutg->SetPoint(1,0.0716332,1.78571);
   cutg->SetPoint(2,0.787966,0.735294);
   cutg->SetPoint(3,2.43553,0.210084);
   cutg->SetPoint(4,4.72779,0.630252);
   cutg->SetPoint(5,6.51862,2.31092);
   cutg->SetPoint(6,8.09456,5.88235);
   cutg->SetPoint(7,7.16332,8.29832);
   cutg->SetPoint(8,6.23209,8.40336);
   cutg->SetPoint(9,4.15473,8.5084);
   cutg->SetPoint(10,2.86533,7.77311);
   cutg->SetPoint(11,0.286533,3.78151);
   cutg->SetPoint(12,0.286533,3.78151);
   cutg->SetPoint(13,0.286533,3.88655);

    //disegno sullo stesso canvas gli eventi, divisi in modo opportuno rispetto al taglio
    TCanvas *c2 = new TCanvas();

    eve->Draw("y:x >> Sig2D","s>0 && !CUTG", "");
    eve->Draw("y:x >> Eff2D", "s>0 && CUTG", "SAME");
    eve->Draw("y:x >> Bkg2D", "s<1 && CUTG", "SAME");
    eve->Draw("y:x >> Pur2D", "s<1 && !CUTG", "SAME");
    
    //calcolo purezza ed efficienza del taglio

    double purezza, efficienza, SS, SB, BS, BB;

    SS = Sig2D->GetEntries();
    SB = Eff2D->GetEntries();
    BS = Pur2D->GetEntries();
    BB = Bkg2D->GetEntries();

    purezza = 100*(SS)/(SS+BS);
    efficienza = 100*(SS)/(SS+SB);

    cout <<"Purezza del taglio: "<<purezza<<"%"<<endl;
    cout <<"Efficienza del taglio: "<<efficienza<<"%"<<endl;
 
    return;
}