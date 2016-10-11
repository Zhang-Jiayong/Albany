#include "CTM_SolutionInfo.hpp"
#include <Albany_AbstractDiscretization.hpp>
#include <PCU.h>

namespace CTM {

    SolutionInfo::SolutionInfo() {
    }

    // get solution vectors
    Teuchos::RCP<Tpetra_MultiVector> SolutionInfo::getOwnedMV() {
        return owned_x;
    }
    
    Teuchos::RCP<Tpetra_MultiVector> SolutionInfo::getGhostMV() {
        return ghost_x;
    }

    // get exporter
    Teuchos::RCP<Tpetra_Export> SolutionInfo::getExporter() {
        return exporter;
    }
    // get importer
    Teuchos::RCP<Tpetra_Import> SolutionInfo::getImporter() {
        return importer;
    }

    // get residual vectors
    Teuchos::RCP<Tpetra_Vector> SolutionInfo::getOwnedResidual() {
        return owned_f;
    }
    //
    Teuchos::RCP<Tpetra_Vector> SolutionInfo::getGhostResidual() {
        return ghost_f;
    }

    // get Jacobian
    Teuchos::RCP<Tpetra_CrsMatrix> SolutionInfo::getOwnedJacobian() {
        return owned_J;
    }
    //
    Teuchos::RCP<Tpetra_CrsMatrix> SolutionInfo::getGhostJacobian() {
        return ghost_J;
    }
    
    void SolutionInfo::gather_x() {
        owned_x->doExport(*ghost_x, *exporter, Tpetra::INSERT);
    }

    void SolutionInfo::scatter_x() {
        ghost_x->doImport(*owned_x, *importer, Tpetra::INSERT);
    }

    void SolutionInfo::gather_f() {
        owned_f->doExport(*ghost_f, *exporter, Tpetra::ADD);
    }

    void SolutionInfo::scatter_f() {
        ghost_f->doImport(*owned_f, *importer, Tpetra::INSERT);
    }

    void SolutionInfo::gather_J() {
        owned_J->doExport(*ghost_J, *exporter, Tpetra::ADD);
    }

    void SolutionInfo::scatter_J() {
        ghost_J->doImport(*owned_J, *importer, Tpetra::INSERT);
    }

    void SolutionInfo::resize(RCP<Albany::AbstractDiscretization> d, bool have_x_dot) {
        auto t0 = PCU_Time();
        int number_vectors = 1;
        if (have_x_dot) number_vectors = 2;
        auto map = d->getMapT();
        auto ghost_map = d->getOverlapMapT();
        auto graph = d->getJacobianGraphT();
        auto ghost_graph = d->getOverlapJacobianGraphT();
        exporter = rcp(new Tpetra_Export(ghost_map, map));
        importer = rcp(new Tpetra_Import(map, ghost_map));
        owned_x = rcp(new Tpetra_MultiVector(map, number_vectors));
        ghost_x = rcp(new Tpetra_MultiVector(ghost_map, number_vectors));
        owned_f = rcp(new Tpetra_Vector(map));
        ghost_f = rcp(new Tpetra_Vector(ghost_map));
        owned_J = rcp(new Tpetra_CrsMatrix(graph));
        ghost_J = rcp(new Tpetra_CrsMatrix(ghost_graph));
        auto t1 = PCU_Time();
        if (!PCU_Comm_Self())
            printf("Solution containers resized in %f seconds", t1 - t0);
    }

}