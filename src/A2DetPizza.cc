// Pizza detector
// Author: Dominik Werthmueller, 2016

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"
#include "G4TessellatedSolid.hh"
#include "CLHEP/Units/SystemOfUnits.h"

#include "A2DetPizza.hh"
#include "A2SD.hh"
#include "A2VisSD.hh"

using namespace CLHEP;

//______________________________________________________________________________
A2DetPizza::A2DetPizza()
{
    // Constructor.

    // init members
    fZPos = 138*cm;
    fIsCheckOverlap = true;
    fRegionPizza = new G4Region("Pizza");
    fPizzaSD = 0;
    fPizzaVisSD = 0;
}

//______________________________________________________________________________
A2DetPizza::~A2DetPizza()
{
    // Destructor.

    if (fRegionPizza) delete fRegionPizza;
    if (fPizzaSD) delete fPizzaSD;
    if (fPizzaVisSD) delete fPizzaVisSD;
}

//______________________________________________________________________________
G4VPhysicalVolume* A2DetPizza::Construct(G4LogicalVolume* motherLogic)
{
    // Construct the detector.

    // global geometry data
    const G4int nPizza = 24;
    G4double phi0 = 360.*deg / nPizza / 2.;
    G4double dphi = 360.*deg / nPizza;
    const G4double scint_thick = 10*mm;

    // set the mother volume
    fMotherLogic = motherLogic;

    // construct the air box
    G4Box* airBox = new G4Box("AirBox", 130.0*cm, 130.0*cm, 4.0*cm);
    fMyLogic = new G4LogicalVolume(airBox, fNistManager->FindOrBuildMaterial("G4_AIR"), "AirBox");
    fMyLogic->SetVisAttributes(G4VisAttributes::Invisible);
    fMyPhysi = new G4PVPlacement(0, G4ThreeVector(0.0*cm, 0.0*cm, fZPos+0.5*scint_thick),
                                                  fMyLogic, "AirBox",
                                                  fMotherLogic, false, 0);
    if (fIsCheckOverlap) CheckOverlapAndAbort(fMyPhysi, "A2DetPizza::Construct()");

    //
    // build the scintillator
    //

    // geometry data (anti-clockwise ordering of surface points)
    const G4int scint_n = 4;
    const G4double scint_x[scint_n] = { 37.734*mm, 727.947*mm, 727.947*mm, 37.734*mm };
    const G4double scint_y[scint_n] = { -4.584*mm, -95.601*mm, 95.601*mm, 4.584*mm };

    // create solid and logical volume
    G4TessellatedSolid* scint = BuildPlanarTessSolid(scint_n, scint_x, scint_y,
                                                     scint_thick, "pizza_scint");
    G4LogicalVolume* scint_log = new G4LogicalVolume(scint,
                                                     fNistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"), // ?
                                                     "pizza_scint");
    scint_log->SetVisAttributes(G4Colour(1, 1, 0));

    // create the physical volumes
    for (G4int i = 0; i < nPizza; i++)
    {
        G4RotationMatrix* rot = new G4RotationMatrix();
        rot->rotateZ(phi0 + i*dphi);
        G4VPhysicalVolume* v = new G4PVPlacement(rot, G4ThreeVector(0.0*cm, 0.0*cm, 0),
                                                 scint_log, "pizza_scint",
                                                 fMyLogic, false, i);
        if (fIsCheckOverlap) CheckOverlapAndAbort(v, "A2DetPizza::Construct()");
    }

    // set-up sensitive detectors
    G4SDManager* sdMan = G4SDManager::GetSDMpointer();
    if (fIsInteractive)
    {
        if (!fPizzaVisSD) fPizzaVisSD = new A2VisSD("PizzaVisSD", nPizza);
        sdMan->AddNewDetector(fPizzaVisSD);
        scint_log->SetSensitiveDetector(fPizzaVisSD);
        fRegionPizza->AddRootLogicalVolume(scint_log);
    }
    else
    {
        if (!fPizzaSD) fPizzaSD = new A2SD("PizzaSD", nPizza);
        sdMan->AddNewDetector(fPizzaSD);
        scint_log->SetSensitiveDetector(fPizzaSD);
        fRegionPizza->AddRootLogicalVolume(scint_log);
    }

    //
    // build the lightguide
    //

    // geometry data (anti-clockwise ordering of surface points)
    const G4int light_guide_n = 10;
    const G4double light_guide_x[light_guide_n] = { 727.947*mm, 837.948*mm, 848.540*mm, 912.113*mm, 940.146*mm,
                                                    940.146*mm, 912.113*mm, 848.540*mm, 837.948*mm, 727.947*mm };
    const G4double light_guide_y[light_guide_n] = { -95.601*mm, -110.201*mm, -110.201*mm, -23.493*mm, -14.905*mm,
                                                    14.905*mm, 23.493*mm, 110.201*mm, 110.201*mm, 95.601*mm };
    const G4double light_guide_thick = 30*mm;

    // create solid and logical volume
    G4TessellatedSolid* light_guide = BuildPlanarTessSolid(light_guide_n, light_guide_x, light_guide_y,
                                                           light_guide_thick, "pizza_light_guide");
    G4LogicalVolume* light_guide_log = new G4LogicalVolume(light_guide,
                                                           fNistManager->FindOrBuildMaterial("A2_PLASTIC"), // ?
                                                           "pizza_light_guide");
    light_guide_log->SetVisAttributes(G4Colour(1, 0, 0));

    // create the physical volumes
    for (G4int i = 0; i < nPizza; i++)
    {
        G4RotationMatrix* rot = new G4RotationMatrix();
        rot->rotateZ(phi0 + i*dphi);
        G4VPhysicalVolume* v = new G4PVPlacement(rot, G4ThreeVector(0.0*cm, 0.0*cm, 0),
                                                 light_guide_log, "pizza_light_guide",
                                                 fMyLogic, false, i);
        if (fIsCheckOverlap) CheckOverlapAndAbort(v, "A2DetPizza::Construct()");
    }

    //
    // build the metal shoes
    //

    // geometry data (anti-clockwise ordering of surface points)
    const G4int shoe_n = 8;
    const G4double shoe_x[shoe_n] = { 728.021*mm, 837.948*mm, 848.540*mm, 918.330*mm,
                                      918.330*mm, 848.540*mm, 837.948*mm, 728.021*mm };
    const G4double shoe_y[shoe_n] = { -95.601*mm, -110.201*mm, -110.201*mm, -109.915*mm,
                                      109.915*mm, 110.201*mm, 110.201*mm, 95.601*mm };
    const G4double shoe_thick = 20*mm;
    const G4double shoe_z = 0.5*light_guide_thick + 0.5*shoe_thick;

    // create solid and logical volume
    G4TessellatedSolid* shoe = BuildPlanarTessSolid(shoe_n, shoe_x, shoe_y,
                                                    shoe_thick, "pizza_shoe");
    G4LogicalVolume* shoe_log = new G4LogicalVolume(shoe,
                                                    fNistManager->FindOrBuildMaterial("G4_Al"), // ?
                                                    "pizza_shoe");
    shoe_log->SetVisAttributes(G4Colour(0, 0, 1));

    // create the physical volumes
    for (G4int i = 0; i < nPizza; i++)
    {
        G4RotationMatrix* rot = new G4RotationMatrix();
        rot->rotateZ(phi0 + i*dphi);
        G4VPhysicalVolume* v = new G4PVPlacement(rot, G4ThreeVector(0.0*cm, 0.0*cm, shoe_z),
                                                 shoe_log, "pizza_shoe_back",
                                                 fMyLogic, false, i);
        if (fIsCheckOverlap) CheckOverlapAndAbort(v, "A2DetPizza::Construct()");
        v = new G4PVPlacement(rot, G4ThreeVector(0.0*cm, 0.0*cm, -shoe_z),
                              shoe_log, "pizza_shoe",
                              fMyLogic, false, i);
        if (fIsCheckOverlap) CheckOverlapAndAbort(v, "A2DetPizza::Construct()");
    }

    //
    // build the PMs
    //

    // geometry data
    G4double innerRadius = 15.5*mm; // ?
    G4double outerRadius = 17.5*mm;
    G4double hz = 0.5*187*mm;
    G4double startAngle = 0*deg;
    G4double spanningAngle = 360*deg;
    G4double light_guide_pos = 940.146*mm;
    G4double center_shift = light_guide_pos+hz;

    // create solid and logical volume
    G4Tubs* pm = new G4Tubs("pizza_pm", innerRadius, outerRadius, hz, startAngle, spanningAngle);
    G4LogicalVolume* pm_log = new G4LogicalVolume(pm,
                                  fNistManager->FindOrBuildMaterial("A2_MUMETAL"), // ?
                                  "pizza_pm");
    pm_log->SetVisAttributes(G4Colour(0, 1, 1));

    // create the physical volumes
    for (G4int i = 0; i < nPizza; i++)
    {
        G4RotationMatrix* rot = new G4RotationMatrix();
        G4double rot_z = phi0 + i*dphi;
        rot->rotateZ(rot_z);
        rot->rotateY(90*deg);
        G4VPhysicalVolume* v = new G4PVPlacement(rot, G4ThreeVector(center_shift*cos(rot_z), -center_shift*sin(rot_z), 0*cm),
                                                 pm_log, "pizza_pm",
                                                 fMyLogic, false, i);
        if (fIsCheckOverlap) CheckOverlapAndAbort(v, "A2DetPizza::Construct()");
    }

    //
    // build the PM protection
    //

    // geometry data
    innerRadius = 25*mm; // ?
    outerRadius = 26*mm; // ? was 22.5 (overlap)
    hz = 0.5*295*mm;
    startAngle = 0*deg;
    spanningAngle = 360*deg;
    G4double shoe_pos = 918.330*mm + 8*mm; // prevents overlap
    center_shift = shoe_pos+hz;

    // create solid and logical volume
    G4Tubs* pm_prot = new G4Tubs("pizza_pm_prot", innerRadius, outerRadius, hz, startAngle, spanningAngle);
    G4LogicalVolume* pm_prot_log = new G4LogicalVolume(pm_prot,
                                                       fNistManager->FindOrBuildMaterial("A2_MUMETAL"), // ?
                                                       "pizza_pm_prot");
    pm_prot_log->SetVisAttributes(G4Colour(1, 0, 1));

    // create the physical volumes
    for (G4int i = 0; i < nPizza; i++)
    {
        G4RotationMatrix* rot = new G4RotationMatrix();
        G4double rot_z = phi0 + i*dphi;
        rot->rotateZ(rot_z);
        rot->rotateY(90*deg);
        G4VPhysicalVolume* v = new G4PVPlacement(rot, G4ThreeVector(center_shift*cos(rot_z), -center_shift*sin(rot_z), 0*cm),
                                                 pm_prot_log, "pizza_pm_prot",
                                                 fMyLogic, false, i);
        if (fIsCheckOverlap) CheckOverlapAndAbort(v, "A2DetPizza::Construct()");
    }

    return fMyPhysi;
}

//______________________________________________________________________________
G4TriangularFacet** A2DetPizza::BuildTriSurface(G4int n, const G4double* x,
                                                const G4double* y, G4double z,
                                                G4bool reverse)
{
    // Build a surface defined by the returned array of 'n' triangular facets
    // using 'n' points with coordinates 'x', 'y' and 'z'.
    // If 'reverse' is true, reverse the order of the points inverting the
    // direction of the normal vector of the surface.

    // check number of points
    if (n < 3)
    {
        G4cout << "A2DetPizza::BuildTriSurface(): At least 3 points needed!" << G4endl;
        return 0;
    }

    // calculate middle point
    G4double mid_x = 0;
    G4double mid_y = 0;
    for (G4int i = 0; i < n; i++)
    {
        mid_x += x[i];
        mid_y += y[i];
    }
    mid_x /= (G4double)n;
    mid_y /= (G4double)n;

    // create triangular facets
    G4TriangularFacet** out = new G4TriangularFacet*[n];
    for (G4int i = 0; i < n; i++)
    {
        // set points
        G4int p1 = i;
        G4int p2 = i == n-1 ? 0 : i+1;

        // create facet
        if (reverse)
        {
            out[i] = new G4TriangularFacet(G4ThreeVector(mid_x, mid_y, z),
                                           G4ThreeVector(x[p2], y[p2], z),
                                           G4ThreeVector(x[p1], y[p1], z),
                                           ABSOLUTE);
        }
        else
        {
            out[i] = new G4TriangularFacet(G4ThreeVector(mid_x, mid_y, z),
                                           G4ThreeVector(x[p1], y[p1], z),
                                           G4ThreeVector(x[p2], y[p2], z),
                                           ABSOLUTE);
        }
    }

    return out;
}

//______________________________________________________________________________
G4QuadrangularFacet** A2DetPizza::BuildDepthFacets(G4int n, const G4double* x,
                                                   const G4double* y, G4double thickness)
{
    // Build an array of of 'n' quadrangular facets using 'n' points with coordinates
    // 'x' and 'y' and use 'thickness' for calculating the z-coordinates.

    // check number of points
    if (n < 2)
    {
        G4cout << "A2DetPizza::BuildDepthFacets(): At least 2 points needed!" << G4endl;
        return 0;
    }

    // create quadrangular facets
    G4QuadrangularFacet** out = new G4QuadrangularFacet*[n];
    for (G4int i = 0; i < n; i++)
    {
        // set points
        G4int p1 = i == n-1 ? 0 : i+1;
        G4int p2 = i;

        // create facet
        out[i] = new G4QuadrangularFacet(G4ThreeVector(x[p1], y[p1], 0.5*thickness),
                                         G4ThreeVector(x[p2], y[p2], 0.5*thickness),
                                         G4ThreeVector(x[p2], y[p2], -0.5*thickness),
                                         G4ThreeVector(x[p1], y[p1], -0.5*thickness),
                                         ABSOLUTE);
    }

    return out;
}

//______________________________________________________________________________
G4TessellatedSolid* A2DetPizza::BuildPlanarTessSolid(G4int n, const G4double* x,
                                                     const G4double* y, G4double thickness,
                                                     const G4String& name)
{
    // Build a planar solid named 'name' that is defined by 'n' points with
    // coordinates 'x' and 'y' and thickness 'thickness'.

    // build the front and back surfaces
    G4TriangularFacet** front_facets = BuildTriSurface(n, x, y, 0.5*thickness);
    G4TriangularFacet** back_facets = BuildTriSurface(n, x, y, -0.5*thickness, true);

    // build the depth facets
    G4QuadrangularFacet** depth_facets = BuildDepthFacets(n, x, y, thickness);

    // create the solid
    G4TessellatedSolid* out = new G4TessellatedSolid(name);

    // add all facets
    for (G4int i = 0; i < n; i++)
    {
        out->AddFacet((G4VFacet*) front_facets[i]);
        out->AddFacet((G4VFacet*) depth_facets[i]);
        out->AddFacet((G4VFacet*) back_facets[i]);
    }

    // mark solid as closed
    out->SetSolidClosed(true);

    return out;
}

//______________________________________________________________________________
void A2DetPizza::CheckOverlapAndAbort(G4VPhysicalVolume* vol, const G4String& location)
{
    // Check if the volume 'vol' has overlaps and abort if that is the case.
    // Use 'location' in the error message to indicate the origin of the
    // problem.

    // check for overlaps
    if (vol->CheckOverlaps())
    {
        G4cout << location << ": Overlap in volume " << vol->GetName() << " detected!" << G4endl;
        exit(1);
    }
}

