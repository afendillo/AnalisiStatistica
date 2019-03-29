#define NRand 1e6
#define NInt 2500
#define bin 200
#define xmax 5
#define xmin 0
#define Range xmax-xmin


using namespace std;




//Poly function
double poly (double x) {

	double y = x*x*x*(x*(3*x + 2) + 3);
	return y;
}

//Trial distribution
double trial (double x) {

	double y = x*x*x*4;
	return y;
}

//Cumulante
double cumul (double x) {

	double y = x*x*x*x;
	return y;
}

//Cumulante invertita
double lumuc (double x) {

	double y = pow(x, 1./4.);
	return y;
}

//Constant PDFunction a-b interval
double cFunc(double x) {
		return ROOT::Math::uniform_pdf  (x,0,1,0);
							
}	
//integrale montecarlo distribuzione uniforme
double integral (double x) {

	double I=0;
	for(int i = 0 ; i<NRand; i++)
		{
			I = I + poly(x);
		}
	return (Range)*I/NRand;
}
//integrale montecarlo distribuzione alternativa
double integral2 (double x) {

	double I=0;
	for(int i = 0 ; i<NRand; i++)
		{
			I = I + poly(x)/trial(x);
		}
	return cumul(xmax)*I/NRand;
}
void Integrale(char stamp='Q')
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
	gStyle->SetOptFit(1111);

//metodo 1

	TRandom3 *RandGenerator = new TRandom3(time(0));
	double Result[NInt];
	double Rmin = 0, Rmax = 0;

	
	for (int i = 0; i<NInt; i++)
	{
		Result[i] = integral(RandGenerator->Rndm()*(Range)+xmin);
		//cout<<Result[i]<<endl;

		if(Result[i]>Rmax) {
			Rmax=Result[i];}

		if(Rmin==0) {
			Rmin=Result[i];}
		else if(Rmin>Result[i]){
			Rmin=Result[i];}
	}

	TH1D *HistoInt = new TH1D ("HistoInt" , "Uniforme" , bin ,Rmin, Rmax);
	HistoInt->GetXaxis()->SetTitle("Area");
	HistoInt->GetYaxis()->SetTitle("Conteggi");
	HistoInt->SetFillColorAlpha(kGreen, 0.30);

	for (int i = 0; i<NInt; i++)
	{
		HistoInt->Fill(Result[i]);
	}


	TCanvas *c1 = new TCanvas();
	c1->SetGrid();
	HistoInt->Draw();

//metodo 2

	double Result2[NInt];
	double ymin = cumul(xmin), ymax = cumul(xmax);
	Rmin = 0;
	Rmax = 0;
	
	for (int i = 0; i<NInt; i++)
	{
		Result2[i] = integral2(lumuc(RandGenerator->Rndm()*(ymax-ymin)+ymin));
		//cout<<Result[i]<<endl;

		if(Result2[i]>Rmax) {
			Rmax=Result2[i];}

		if(Rmin==0) {
			Rmin=Result2[i];}
		else if(Rmin>Result2[i]){
			Rmin=Result2[i];}
	}

	TH1D *HistoInt2 = new TH1D ("HistoInt2" , "Cubica" , bin ,Rmin, Rmax);
	HistoInt2->GetXaxis()->SetTitle("Area");
	HistoInt2->GetYaxis()->SetTitle("Conteggi");
	HistoInt2->SetFillColorAlpha(kGreen, 0.30);

	for (int i = 0; i<NInt; i++)
	{
		HistoInt2->Fill(Result2[i]);
	}


	TCanvas *c2 = new TCanvas();
	c2->SetGrid();
	HistoInt2->Draw();
	

	return;	
}

