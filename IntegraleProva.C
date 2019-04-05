#define NRand 1e3
#define NInt 2.5e3
#define bin 50
#define xmax 5
#define xmin 0
#define Range xmax-xmin


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
//Poly function
double poly (double x) {
	return x*x*x*(x*(3*x + 2) + 3);
}

//Trial distribution
double trial (double x) {
	return  x*x*x*x*x*6;
}

//Cumulante
double cumul (double x) {
	return x*x*x*x*x*x;
}

//Cumulante invertita
double lumuc (double x) {
	return pow(x, 1./6.);
}

//Trial distribution
double trials (double x) {
	return  x*x*x*4;
}

//Cumulante
double cumuls (double x) {
	return x*x*x*x;
}

//Cumulante invertita
double slumuc (double x) {
	return pow(x, 1./4.);
}
//Constant PDFunction a-b interval
double cFunc(double x) {
		return ROOT::Math::uniform_pdf  (x,0,1,0);						
}	
//integrale montecarlo distribuzione uniforme
// double integral (double* rand) {

// 	double I=0, x;
// 	for(int i = 0 ; i<NRand; i++)
// 		{
// 			x=rand[i]*(Range)+xmin;
// 			I = I + poly(x);
// 		}
// 	return (Range)*I/NRand;
// }
//integrale montecarlo distribuzione alternativa
// double integral2 (double* rand) {

// 	double ymin = cumul(xmin), ymax = cumul(xmax);
// 	double I=0, x;
// 	for(int i = 0 ; i<NRand; i++)
// 		{
// 			x=lumuc(rand[i]*(ymax-ymin)+ymin);
// 			I = I + poly(x)/trial(x);
// 		}
// 	return cumul(xmax)*I/NRand;
// }
void Integrale()
{	
	Reset();
	gStyle->SetOptFit(1111);
//metodo 1

	TRandom3 *RandGenerator = new TRandom3(time(0));
	double Result , Rmin = 0, Rmax = 0 , I ,J, IS, Is, x , y , xS;
	double Result2, ymin = cumul(xmin), ymax = cumul(xmax), ymid = cumul(xmax/2.);
	double ymids = cumuls(xmax/2), ymins = cumuls(xmin), ymaxs = cumuls(xmax);

	TH1D *HistoInt = new TH1D ("HistoInt" , "Uniforme" , bin ,Rmin, Rmax);
	HistoInt->GetXaxis()->SetTitle("Area");
	HistoInt->GetYaxis()->SetTitle("Conteggi");
	HistoInt->SetFillColorAlpha(kGreen, 0.30);

	TH1D *HistoInt2 = new TH1D ("HistoInt2" , "Cubica" , bin ,Rmin, Rmax);
	HistoInt2->GetXaxis()->SetTitle("Area");
	HistoInt2->GetYaxis()->SetTitle("Conteggi");
	HistoInt2->SetFillColorAlpha(kGreen, 0.30);

	TH1D *HistoInt3 = new TH1D ("HistoInt3" , "Samplig Stratificato" , bin ,Rmin, Rmax);
	HistoInt3->GetXaxis()->SetTitle("Area");
	HistoInt3->GetYaxis()->SetTitle("Conteggi");
	HistoInt3->SetFillColorAlpha(kGreen, 0.30);

	TF1 *f = new TF1("f","gaus(0)" , 0 , 10);
	f->SetParameter(0 , 1);
	f->SetParameter(1 , 9531.25);
	f->SetParameter(2 , 1);

	double RangeS=Range/2, xminS=0;
	cout<<"\nStarting....\n";
	for (int i = 0; i<NInt; i++)
	{
		I=0;J=0;IS=0;
		if(i%(int)(NInt/100)==0) 
		{
			for (int k = 0 ; k<i/(int)(NInt/100)+1; k++) cout << "*" ;
			cout<<" "<< i/(int)(NInt/100)+1<< "% \r";
			cout<<flush;
		}

		if(i<NInt/2) xminS=xmin;
		else xminS=xmax/2;

		for(int j = 0 ; j<NRand; j++)
		{
			
			//x=(RandGenerator->Rndm())*(Range)+xmin;
			//I = I + poly(x);

			y=lumuc((RandGenerator->Rndm())*(ymax-ymin)+ymin);
			J = J + poly(y)/trial(y);
			
			if(i<NInt/2){
				xS=slumuc((RandGenerator->Rndm())*(ymids-ymins)+ymins);
				Is = Is + poly(xS)/trials(xS);
			}
			else{
				xS=lumuc((RandGenerator->Rndm())*(ymax-ymid)+ymid);
				IS = IS + poly(xS)/trial(xS);
			}
			
		}
		HistoInt2->Fill( cumul(xmax)*J/NRand);
		//HistoInt->Fill((Range)*I/NRand);
		HistoInt3->Fill(((cumul(xmax)-cumul(xmax/2))*IS/(NRand/2))+((cumuls(xmax/2)-cumuls(xmin))*Is/(NRand/2)));
	}
	cout<<endl;

/*	TCanvas *c1 = new TCanvas();
	c1->SetGrid();
	HistoInt->Draw();
	HistoInt->Fit("f", "Q");
*/
	TCanvas *c2 = new TCanvas();
	c2->SetGrid();
	HistoInt2->Draw();
	HistoInt2->Fit("f" , "Q");

	TCanvas *c3 = new TCanvas();
	c3->SetGrid();
	HistoInt3->Draw();
	HistoInt3->Fit("f" , "Q");

	return;	
}

