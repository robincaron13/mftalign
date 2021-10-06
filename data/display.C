{

  //TGeoManager::Import("geofile_mft_orig.root");
  TGeoManager::Import("o2sim_geometry_misaligned.root");
    
  //TGeoVolume *mft = gGeoManager->GetVolume("MFT_0");

  TGeoVolume *mftgeo = gGeoManager->GetTopVolume();
  TGeoVolume *mft2 = gGeoManager->GetVolume("MFT_H_1");

  TGeoVolume *disk0 = gGeoManager->GetVolume("MFT_D_0_0");
  TGeoVolume *disk1 = gGeoManager->GetVolume("MFT_D_1_1");
  TGeoVolume *disk3 = gGeoManager->GetVolume("MFT_D_3_3");

    
  //gGeoManager->SetVisLevel(10);
  mftgeo->Draw("ogl");
  //mft2->Draw("ogl same");
  //disk0->Draw("ogl");
  //disk1->Draw("ogl");
  //disk3->Draw("ogl same");

    
}

