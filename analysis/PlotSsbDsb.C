/*
Copyright 2026 Shun Fukagawa, Tsukasa Aso

  Licensed under the Apache License, Version 2.0.
  You may not use this file except in compliance with the License.
  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0
*/


#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TStyle.h"

#include "TApplication.h"
#include "TROOT.h"

#include <iostream>
#include <string>

void PlotSsbDsb(TString filename){

    gStyle->SetOptStat(0); // disable statistics box

    // Create canvas
    TCanvas *canvas = new TCanvas("canvas", "Graph", 800, 800);
    canvas->SetLogx();  // Set x-axis to logarithmic scale
    canvas->SetLogy();  // Set y-axis to logarithmic scale
    //canvas->SetGrid();

    auto gr0 = new TGraphErrors(filename, "%lg %lg %lg %lg", " ");
    gr0->SetMarkerStyle(8);               // Circle
    gr0->SetMarkerColor(kBlue-3);         // Blue
    gr0->SetMarkerSize(1.5);
    gr0->SetLineWidth(2);
    gr0->SetLineColor(kBlack);            // Set line color to black
    gr0->SetTitle("");

    // Set graph drawing range
    gr0->GetXaxis()->SetLimits(1e-1, 1e3);
    // gr0->GetYaxis()->SetLimits(1e0 , 1e2);
    gr0->GetYaxis()->SetRangeUser(1e0, 1e2);
    gr0->GetXaxis()->SetNoExponent(kFALSE);    // Allow exponential notation
    gr0->GetYaxis()->SetNoExponent(kFALSE);

    // Set axis labels
    gr0->GetXaxis()->SetTitleFont(62);
    gr0->GetYaxis()->SetTitleFont(62);
    gr0->GetXaxis()->SetTitleOffset(1.3);
    gr0->GetXaxis()->CenterTitle();
    gr0->GetYaxis()->CenterTitle();
    gr0->GetXaxis()->SetTitle("LET (keV/#mum)");
    gr0->GetYaxis()->SetTitle("SSB/DSB");

    gPad->SetTickx(1);   // Draw ticks on both top and bottom of X axis
    gPad->SetTicky(1);   // Draw ticks on both left and right of Y axis

    // Draw graph
    gr0->Draw("AP");
    canvas->Update();

    // Display axis tick labels in exponential notation
    TAxis *xa = gr0->GetXaxis();
    xa->ChangeLabel(1, -1, -1, -1, -1, -1, "10^{-1}");
    xa->ChangeLabel(2, -1, -1, -1, -1, -1, "10^{0}");
    xa->ChangeLabel(3, -1, -1, -1, -1, -1, "10^{1}");
    xa->ChangeLabel(4, -1, -1, -1, -1, -1, "10^{2}");
    xa->ChangeLabel(5, -1, -1, -1, -1, -1, "10^{3}");

    TAxis *ya = gr0->GetYaxis();
    ya->ChangeLabel(1, -1, -1, -1, -1, -1, "10^{0}");
    ya->ChangeLabel(2, -1, -1, -1, -1, -1, "10^{1}");
    ya->ChangeLabel(3, -1, -1, -1, -1, -1, "10^{2}");

    canvas->Modified();
    canvas->Update();
    
    // display legend
    //canvas->BuildLegend();
    auto leg1 = new TLegend(0.52, 0.61, 0.89, 0.89);

    leg1->SetBorderSize(0);
    leg1->Draw();

    // Update canvas
    canvas->Update();
    std::string res;
    std::cout << "Do you want to save the image? [y/n]: ";
    std::cin >> res;
    if (res == "y") {
        TString outFilename;
        std::cout << "Enter the output file name: ";
        std::cin >> outFilename;
        canvas->SaveAs(outFilename.Data());
    }
}
