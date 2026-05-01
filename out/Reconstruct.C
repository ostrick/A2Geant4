
//Global definitions
Double_t vdrift = 354; //cm/ms: 1/10 of what is in A2DriftModel::SetConstants
Double_t qeslope = -0.000208; //from output of Calibrate.C
Double_t qeintercept = -0.89; //from output of Calibrate.C

//Attach an anode section ID number to an x coordinate for theta reconstruction
Double_t GetX(Int_t secID){
	//based on anode numbering scheme in g4
        Int_t row1[16] = {1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61};
        Int_t row2[16] = {2,6,10,14,18,22,26,30,34,38,42,46,50,54,58,62};
        Int_t row3[16] = {3,7,11,15,19,23,27,31,35,39,43,47,51,55,59,63};
        Int_t row4[16] = {4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64};
        for (int i=0;i<16;i++){ //over all angular sections in one ring
		//these numbers need to be changed if testing different anode designs!
                if (secID==row1[i])return 10; //assign the radius of that ring
                if(secID==row2[i])return 20;
                if(secID==row3[i])return 30;
                if(secID==row4[i])return 40;
        }
        if (secID==65)return 5; //central ring
        return 0; //central pad
}

//Reconstruct recoil polar angle
Double_t GetTheta(Double_t mintime, Double_t maxtime, Double_t minx, Double_t maxx){
	Double_t deltax, deltat, theta;
	theta=0;
	deltax=maxx-minx;
	deltat=maxtime-mintime;
	if(deltax!=0 && deltat!=0)theta=90-atan((deltat*vdrift*10)/deltax)*(180/3.14);
	return theta;
}

//Reconstruct event position
Double_t GetZ(Double_t time){
	Double_t z_pos =-11.5 +vdrift *time;
	return z_pos;
}

//Reconstruct recoil energy
Double_t GetEnergy(Double_t charge){
	Double_t energy = qeintercept + (qeslope)*charge;
	return energy;
}

//Draw the reconstruction output
void Draw(TFile* output){
	//draws fractional error distributions
	//and reconstructed vs true graphs
	TTree *recon = (TTree*)output->Get("recon");
	TCanvas *c4 = new TCanvas("c4","Reconstruction Output");
	c4->Divide(3,2);
	c4->cd(1);
	recon->Draw("(theta_rec-theta_true)/theta_true","theta_rec!=0");
	TH1F *acc = (TH1F*)gPad->GetPrimitive("htemp");
	acc->SetName("accuracy");
	acc->SetTitle("Fractional Error: Angular Reconstruction");
	acc->GetXaxis()->SetTitle("(Reconstructed-True)/True");
	acc->SetLineColor(kBlue+3);
	acc->SetLineWidth(3);
	acc->SetFillColor(kBlue-7);
	c4->cd(2);
	recon->Draw("(energy_rec-energy_true)/energy_true","energy_rec>0");
	TH1F *acce = (TH1F*)gPad->GetPrimitive("htemp");
	acce->SetName("accuracy");
	acce->SetTitle("Fractional Error: Energy Reconstruction");
	acce->GetXaxis()->SetTitle("(Reconstructed-True)/True");
	acce->SetLineColor(kBlue+3);
	acce->SetLineWidth(3);
	acce->SetFillColor(kBlue-7);
	c4->cd(3);
	recon->Draw("(z_rec-z_true)/(z_true+11.5)","");
	TH1F *accz = (TH1F*)gPad->GetPrimitive("htemp");
	accz->SetName("accuracy");
	accz->SetTitle("Accuracy of Position Reconstruction");
	accz->GetXaxis()->SetTitle("(Reconstructed-True)/True");
	accz->SetLineColor(kBlue+3);
	accz->SetLineWidth(3);
	accz->SetFillColor(kBlue-7);
	c4->cd(4);
	recon->Draw("theta_rec:theta_true","theta_rec!=0 && z_rec!=-11.5 && energy_rec>0","COLZ");
	TH2F *theta = (TH2F*)gPad->GetPrimitive("htemp");
	theta->SetTitle("Recoil Polar Angle (deg)");
	theta->GetXaxis()->SetTitle("True");
	theta->GetYaxis()->SetTitle("Reconstructed");
	c4->cd(5);
	recon->Draw("z_rec:z_true","abs(z_rec)<11.5 && energy_rec>0","");
	TH2F *zco = (TH2F*)gPad->GetPrimitive("htemp");
	zco->SetTitle("Event Z Position (cm)");
	zco->GetXaxis()->SetTitle("True");
	zco->GetYaxis()->SetTitle("Reconstructed");
	c4->cd(6);
	recon->Draw("energy_rec:energy_true","energy_rec>0 && z_rec!=-11.5","");
	TH2F *ke = (TH2F*)gPad->GetPrimitive("htemp");
	ke->SetTitle("Recoil Energy (MeV)");
	ke->GetXaxis()->SetTitle("True");
	ke->GetYaxis()->SetTitle("Reconstructed");
}


//Main function: read data, reconstruct, save to file, and draw
void Reconstruct(TString filename){
	//read in data file and tree
	TFile *in = new TFile(filename);
	TTree *h12 = (TTree*)in->Get("h12");
	//define variables for input tree information
	Int_t ntpc;
	Int_t *itpc = new Int_t[100];
	Float_t *qtpc = new Float_t[100];
	Float_t *ttpc = new Float_t[100];
	Float_t *vertex = new Float_t[3];
	Float_t *klab = new Float_t[3];
	Float_t dircos[100][3];
	//set branch addresses for input tree
	h12->SetBranchAddress("ntpc",&ntpc);
	h12->SetBranchAddress("itpc",itpc);
	h12->SetBranchAddress("qtpc",qtpc);
	h12->SetBranchAddress("ttpc",ttpc);
	h12->SetBranchAddress("klab",klab);
	h12->SetBranchAddress("vertex",vertex);
	h12->SetBranchAddress("dircos",&dircos);
	//define variables to be used in computations
	Int_t first_sec, last_sec;
	Float_t first_time, last_time, charge;
	Int_t max = h12->GetEntries();
	//define variables for output tree information
	Float_t theta_true, z_true, energy_true;
	Float_t theta_rec, z_rec, energy_rec;
	Double_t hits[100][2];
	Float_t deltat;
	//define output file and tree
	TFile *out = new TFile("Reconstructed.root","RECREATE","Reconstructed TPC Data");
	TTree *goat = new TTree("recon","Reconstructed TPC Data");
	//set branch addresses for output tree
	goat->Branch("theta_rec",&theta_rec,"theta_rec/F"); //reconstructed recoil theta
	goat->Branch("theta_true",&theta_true,"theta_true/F"); //simulated recoil theta
	goat->Branch("z_rec",&z_rec,"z_rec/F"); //reconstructed z coordinate
	goat->Branch("z_true",&z_true,"z_true/F"); //simulated z coordinate
	goat->Branch("energy_rec",&energy_rec,"energy_rec/F"); //reconstructed recoil energy
	goat->Branch("energy_true",&energy_true,"energy_true/F"); //simulated recoil energy
	//goat->Branch("dt",&deltat,"dt/F");
	Double_t mintime, maxtime, minx, maxx;
	//read input file, do computations, and write output file
	for (int i=0;i<max;i++){
		h12->GetEntry(i);
		charge=0; //initialize
		//sum charge over all anode sections
		for (int j=0;j<ntpc;j++){
			charge+=qtpc[j];
		}
		//get two data points for theta
		//Note index 0 = most recently saved by simulation = last to occur
		//index ntpc-1 = first saved by simulation = first to occur
		//so point 1 = "min" = tpc[ntpc-1], point 2 = "max" = tpc[0]
		maxtime=ttpc[0];
		maxx=GetX(itpc[0]);
		mintime=ttpc[ntpc-1];
		minx=GetX(itpc[ntpc-1]);
		//try not to use central pad for min if hits in both ring and one other section
			
		if (minx==0 && maxx !=5){
			for (int j=0; j<ntpc; j++){
				if (itpc[j]==65){
					minx=5;
					mintime=ttpc[j];
				}
			}
		}
		
	
		//reconstruct variables
		z_rec=GetZ(ttpc[ntpc-1]); //central section: along beam line
		energy_rec=GetEnergy(charge);
		theta_rec=GetTheta(mintime,maxtime,minx,maxx);	
		//get true variables
		theta_true=acos(dircos[1][2])*180/3.14;
		z_true=vertex[2];
		energy_true=klab[1];
		//write to tree
		goat->Fill();
	}
	out->Write(); //write the goat tree to a file and make the file
	Draw(out);
}
