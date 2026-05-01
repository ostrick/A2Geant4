/*** Calibrate ACTAM TPC 
 * This macro is run using the ouput of the g4 macro TPCcalibrate.mac.
 * It generates the fit parameters of the charge-energy relation
 * for use in event reconstruction.
 * **** ACP 2021 *****/



//This function is called by the main function to find the calibration points
//Reads a g4 output file and returns histogram of total charge deposited in each event
TH1D* GetPoint(TString name){
	TFile *file = new TFile(name);
	TTree *h12 = (TTree*)file->Get("h12");
	TH1D *energy_cal = new TH1D("calibrateE","Energy Calibration",45000,-45000,0); //can fix bins later
	//set some branches for tree reading
	Int_t ntpc;
	Float_t *qtpc = new Float_t[100];
	h12->SetBranchAddress("ntpc",&ntpc);
	h12->SetBranchAddress("qtpc",qtpc);
	//define variables for histogram
	Double_t charge;
	//fill a histogram with the total charge in each event
	Int_t max = h12->GetEntries();
	for (int i=0;i<max;i++){
		h12->GetEntry(i);
		charge=0;
		if (ntpc!=0){
			//sum total charge
			for (int j=0; j<ntpc; j++){
				charge+=qtpc[j];
			}
		energy_cal->Fill(charge);
		}
	}
	//return the histogram
	return energy_cal;
}

//Main function
//Generates fit parameters for a linear relationship
//Recoil KE as a function of charge
void Calibrate(){
	//main function - generates fit parameters
	//define variables necessary for the regression
	Float_t e_mev[16], e_err[16], q_avg[16], q_err[16]; //charge, energy, errors
	//read in the files to assign data
	Int_t e1, e0;
	for(int i=2; i<23; i++){ //23 calibration points
	        e_mev[i-2] = 0.5*i; //every 2 MeV from 2 to 12
       		e_err[i-2] = 0.1; 
		if (i%2==0){
			e0=0;
			e1=i/2;
		}
		else{
		       	e0=5;
			e1=(i-1)/2;
		}
		TString qe_file = Form("Cal_%i_%i.root",e1, e0); //get appropriate file
 		TH1D *qe_point = GetPoint(qe_file); //get correct histogram
		q_avg[i-2]=qe_point->GetMean(); //assign data from histogram
		q_err[i-2]=qe_point->GetRMS();
	}

	//create graph from this data
	TGraphErrors *qegraph = new TGraphErrors(14,q_avg,e_mev,q_err,e_err);
	//optional: draw the graph and fit 
	TCanvas *c0 = new TCanvas("c0","Calibration Data Fits");
	qegraph->SetMarkerColor(1);
        qegraph->SetMarkerSize(1.5);
        qegraph->SetMarkerStyle(21);
        qegraph->SetLineWidth(4);
	qegraph->SetTitle("Energy Calibration");
	qegraph->GetXaxis()->SetTitle("Charge Deposited (e)");
	qegraph->GetYaxis()->SetTitle("Recoil Energy (MeV)");
	qegraph->Draw("AP");
	//create linear fit
	TF1 *qefit = new TF1("qefit","[0]+[1]*x");
	//and fit to data!
	qegraph->Fit(qefit);
}

