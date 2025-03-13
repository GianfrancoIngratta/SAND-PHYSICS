int Get_graphite_polyp_mass_ratio() {
    // Carica il file e crea il gestore TGeoManager
    auto g = TGeoManager::Import("/storage/gpfs_data/neutrino/users/gi/dunendggd/SAND_opt3_DRIFT1.root");
    
    // Verifica se il file è stato caricato correttamente
    if (!g) {
        std::cerr << "Error: Could not load geometry file!" << std::endl;
        return -1;
    }
    
    // Definisci i nomi dei volumi da cercare
    std::vector<std::string> volumes = {"C3H6Target_X0", "C3H6Target_X1", "C3H6Target_A", "C3H6Target_B", "C3H6Target_C"};
    
    // Itera su ciascun nome di volume
    for (const auto& volName : volumes) {
        // Trova il volume nel file di geometria
        TGeoVolume* v = g->FindVolumeFast(volName.c_str());
        
        // Verifica se il volume è stato trovato
        if (!v) {
            std::cerr << "Warning: Volume " << volName << " not found!" << std::endl;
            continue;
        }
        
        // Ottieni la forma e verifica se è un parallelepipedo (TGeoBBox)
        TGeoShape* shape = v->GetShape();
        if (shape->IsA() == TGeoBBox::Class()) {
            TGeoBBox* box = dynamic_cast<TGeoBBox*>(shape);
            if (box) {
                // Stampa le dimensioni del parallelepipedo
                std::cout << "Volume: " << volName << std::endl;
                std::cout << "  Dimensions (half-lengths): " 
                          << "dz = " << box->GetDX() << " cm, "
                          << "dy = " << box->GetDY() << " cm, "
                          << "dx = " << box->GetDZ() << " cm" << std::endl;
                
                // Ottieni il materiale associato e la sua densità
                TGeoMaterial* material = v->GetMedium()->GetMaterial();
                double density = material->GetDensity();  // densità in g/cm^3
                std::cout << "  Density: " << density << " g/cm^3" << std::endl;
            }
        } else {
            std::cerr << "Warning: Volume " << volName << " does not have a TGeoBBox shape." << std::endl;
        }
    }

    return 0;
}