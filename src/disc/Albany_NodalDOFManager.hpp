//*****************************************************************//
//    Albany 3.0:  Copyright 2016 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

#ifndef ALBANY_NODAL_DOF_MANAGER_HPP
#define ALBANY_NODAL_DOF_MANAGER_HPP

#include "Albany_ScalarOrdinalTypes.hpp"
#include "Albany_ThyraTypes.hpp"
#include "Albany_GlobalLocalIndexer.hpp"

#include "Teuchos_RCP.hpp"

namespace Albany {

class NodalDOFManager {
public:
  NodalDOFManager ()  = default;

  void setup (const int numComponents, const LO numLocalNodes,
              const GO maxGlobalNodeID, const bool interleaved) {
    m_numComponents = numComponents;
    m_numLocalNodes = numLocalNodes;
    m_maxGlobalNodeID = maxGlobalNodeID;
    m_interleaved = interleaved;
  }

  void setup (const int numComponents, const Teuchos::RCP<const Thyra_VectorSpace>& node_vs,
              const bool interleaved) {
    const auto indexer = createGlobalLocalIndexer(node_vs);
    const int numLocalNodes = indexer->getNumLocalElements();
    const int maxGlobalGID  = indexer->getMaxGlobalGID();
    setup(numComponents,numLocalNodes,maxGlobalGID,interleaved);
  }

  inline LO getLocalDOF(LO inode, int icomp) const {
    if (m_interleaved) {
      return inode*m_numComponents + icomp;
    } else {
      return inode + m_numLocalNodes*icomp;
    }
  }
  inline GO getGlobalDOF(GO node, int icomp) const {
    if (m_interleaved) {
      return node*m_numComponents + icomp;
    } else {
      return node + m_numLocalNodes*icomp;
    }
  }

  int numComponents() const {
    return m_numComponents;
  }

private:
  int  m_numComponents;
  LO   m_numLocalNodes;
  GO   m_maxGlobalNodeID;
  bool m_interleaved;
};

} // namespace Albany

#endif // ALBANY_NODAL_DOF_MANAGER_HPP
