/// \ingroup base
/// \class ttk::PersistenceDiagram
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \author Julien Tierny <julien.tierny@lip6.fr>
/// \date September 2016.
///
/// \brief TTK processing package for the computation of persistence diagrams.
///
/// This package computes the persistence diagram of the extremum-saddle pairs
/// of an input scalar field. The X-coordinate of each pair corresponds to its
/// birth, while its smallest and highest Y-coordinates correspond to its birth
/// and death respectively.
///
/// In practice, each extremity of a persistence pair is represented by its
/// vertexId and critical type. Based on that, the persistence of the pair
/// and its 2D embedding can easily be obtained.
///
/// Persistence diagrams are useful and stable concise representations of the
/// topological features of a data-set. It is useful to fine-tune persistence
/// thresholds for topological simplification or for fast similarity
/// estimations for instance.
///
/// \b Related \b publication \n
/// "Computational Topology: An Introduction" \n
/// Herbert Edelsbrunner and John Harer \n
/// American Mathematical Society, 2010
///
/// Five backends are available for the computation:
///
///  1) FTM \n
/// \b Related \b publication \n
/// "Task-based Augmented Contour Trees with Fibonacci Heaps"
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny
/// IEEE Transactions on Parallel and Distributed Systems, 2019
///
///  2) Progressive Approach \n
/// \b Related \b publication \n
/// "A Progressive Approach to Scalar Field Topology" \n
/// Jules Vidal, Pierre Guillou, Julien Tierny\n
/// IEEE Transactions on Visualization and Computer Graphics, 2021
///
/// 3) Discrete Morse Sandwich (default) \n
/// \b Related \b publication \n
/// "Discrete Morse Sandwich: Fast Computation of Persistence Diagrams for
/// Scalar Data -- An Algorithm and A Benchmark" \n
/// Pierre Guillou, Jules Vidal, Julien Tierny \n
/// IEEE Transactions on Visualization and Computer Graphics, 2023.\n
/// arXiv:2206.13932, 2023.\n
/// Fast and versatile algorithm for persistence diagram computation.
///
/// 4) Approximate Approach \n
/// \b Related \b publication \n
/// "Fast Approximation of Persistence Diagrams with Guarantees" \n
/// Jules Vidal, Julien Tierny\n
/// IEEE Symposium on Large Data Visualization and Analysis (LDAV), 2021
///
/// 5) Persistent Simplex \n
/// This is a textbook (and very slow) algorithm, described in
/// "Algorithm and Theory of Computation Handbook (Second Edition)
/// - Special Topics and Techniques" by Atallah and Blanton on page 97.
///
/// \sa ttkPersistenceDiagram.cpp %for a usage example.
///
/// \b Online \b examples: \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/1manifoldLearning/">1-Manifold
///   Learning example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/1manifoldLearningCircles/">1-Manifold
///   Learning Circles example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/2manifoldLearning/">
///   2-Manifold Learning example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/BuiltInExample1/">BuiltInExample1
///   </a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/clusteringKelvinHelmholtzInstabilities/">
///   Clustering Kelvin Helmholtz Instabilities example</a> \n
///   - <a href="https://topology-tool-kit.github.io/examples/ctBones/">CT Bones
///   example</a> \n
///   - <a href="https://topology-tool-kit.github.io/examples/dragon/">Dragon
///   example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/harmonicSkeleton/">
///   Harmonic Skeleton example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/imageProcessing/">Image
///   Processing example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/interactionSites/">
///   Interaction sites</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/karhunenLoveDigits64Dimensions//">Karhunen-Love
///   Digits 64-Dimensions example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morsePersistence/">Morse
///   Persistence example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/morseSmaleQuadrangulation/">Morse-Smale
///   Quadrangulation example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 0 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 1 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 2 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 3 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceClustering0/">Persistence
///   clustering 4 example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceDiagramClustering/">Persistence
///   Diagram Clustering example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/persistenceDiagramDistance/">Persistence
///   Diagram Distance example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/tectonicPuzzle/">Tectonic
///   Puzzle example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/tribute/">Tribute
///   example</a> \n
///   - <a
///   href="https://topology-tool-kit.github.io/examples/uncertainStartingVortex/">
///   Uncertain Starting Vortex example</a> \n

#pragma once

// base code includes
#include <ApproximateTopology.h>
#include <DiscreteMorseSandwich.h>
#include <DiscreteMorseSandwichMPI.h>
#include <FTMTreePP.h>
#include <PersistenceDiagramUtils.h>
#include <PersistentSimplexPairs.h>
#include <ProgressiveTopology.h>
#include <Triangulation.h>
#include <psort.h>
#include <string>

namespace ttk {

  namespace persistenceSort {
    inline bool comp(const PersistencePair a, const PersistencePair b) {
      return a.birth.offset < b.birth.offset;
    };

    inline bool oppositeComp(const PersistencePair a, const PersistencePair b) {
      return a.birth.offset > b.birth.offset;
    };
  } // namespace persistenceSort

  /**
   * Compute the persistence diagram of a function on a triangulation.
   * TTK assumes that the input dataset is made of only one connected component.
   */
  class PersistenceDiagram : virtual public Debug {

  public:
    enum class BACKEND {
      FTM = 0,
      PROGRESSIVE_TOPOLOGY = 1,
      DISCRETE_MORSE_SANDWICH = 2,
      APPROXIMATE_TOPOLOGY = 3,
      PERSISTENT_SIMPLEX = 4,
    };

    PersistenceDiagram();

    inline void setBackend(const BACKEND be) {
      this->BackEnd = be;
    }

    inline void setComputeMinSad(const bool data) {
      this->dms_.setComputeMinSad(data);
    }
    inline void setComputeSadSad(const bool data) {
      this->dms_.setComputeSadSad(data);
    }
    inline void setComputeSadMax(const bool data) {
      this->dms_.setComputeSadMax(data);
    }

    /**
     * @brief Complete a ttk::DiagramType instance with scalar field
     * values (useful for persistence) and 3D coordinates of critical vertices
     */
    template <typename scalarType, typename triangulationType>
    void
      augmentPersistenceDiagram(std::vector<PersistencePair> &persistencePairs,
                                const scalarType *const scalars,
                                const triangulationType *triangulation);

    ttk::CriticalType getNodeType(ftm::FTMTree_MT *tree,
                                  ftm::TreeType treeType,
                                  const SimplexId vertexId) const;

    void sortPersistenceDiagram(std::vector<PersistencePair> &diagram,
                                const SimplexId *const offsets) const;

    template <typename scalarType>
    int computeCTPersistenceDiagram(
      ftm::FTMTreePP &tree,
      const std::vector<
        std::tuple<ttk::SimplexId, ttk::SimplexId, scalarType, bool>> &pairs,
      std::vector<PersistencePair> &diagram) const;

    /**
     * @pre For this function to behave correctly in the absence of
     * the VTK wrapper, ttk::preconditionOrderArray() needs to be
     * called to fill the @p inputOffsets buffer prior to any
     * computation (the VTK wrapper already includes a mechanism to
     * automatically generate such a preconditioned buffer).
     * @see examples/c++/main.cpp for an example use.
     */
    template <typename scalarType, class triangulationType>
    int execute(std::vector<PersistencePair> &CTDiagram,
                const scalarType *inputScalars,
                const size_t scalarsMTime,
                const SimplexId *inputOffsets,
                const triangulationType *triangulation);

    template <typename scalarType, class triangulationType>
    int executeFTM(std::vector<PersistencePair> &CTDiagram,
                   const scalarType *inputScalars,
                   const SimplexId *inputOffsets,
                   const triangulationType *triangulation);

    template <class triangulationType>
    int executeProgressiveTopology(std::vector<PersistencePair> &CTDiagram,
                                   const SimplexId *inputOffsets,
                                   const triangulationType *triangulation);
    template <typename scalarType, class triangulationType>
    int executeApproximateTopology(std::vector<PersistencePair> &CTDiagram,
                                   const scalarType *inputScalars,
                                   const triangulationType *triangulation);

    template <class triangulationType>
    int executePersistentSimplex(std::vector<PersistencePair> &CTDiagram,
                                 const SimplexId *inputOffsets,
                                 const triangulationType *triangulation);

    template <typename scalarType, class triangulationType>
    int executeDiscreteMorseSandwich(std::vector<PersistencePair> &CTDiagram,
                                     const scalarType *inputScalars,
                                     const size_t scalarsMTime,
                                     const SimplexId *inputOffsets,
                                     const triangulationType *triangulation);

    template <class triangulationType>
    void checkProgressivityRequirement(const triangulationType *triangulation);

    template <class triangulationType>
    void checkManifold(const triangulationType *const triangulation);

    inline void
      preconditionTriangulation(AbstractTriangulation *triangulation) {
      if(triangulation) {
        triangulation->preconditionBoundaryVertices();
        if(this->BackEnd == BACKEND::FTM
           || this->BackEnd == BACKEND::PROGRESSIVE_TOPOLOGY
           || this->BackEnd == BACKEND::APPROXIMATE_TOPOLOGY) {
          contourTree_.setDebugLevel(debugLevel_);
          contourTree_.setThreadNumber(threadNumber_);
          contourTree_.preconditionTriangulation(triangulation);
        }
        if(this->BackEnd == BACKEND::DISCRETE_MORSE_SANDWICH) {
          dms_.setDebugLevel(debugLevel_);
          dms_.setThreadNumber(threadNumber_);
          dms_.preconditionTriangulation(triangulation);
          triangulation->preconditionManifold();
        }
        if(this->BackEnd == BACKEND::PERSISTENT_SIMPLEX
           || this->BackEnd == BACKEND::DISCRETE_MORSE_SANDWICH) {
          psp_.preconditionTriangulation(triangulation);
        }
      }
    }

    inline void setOutputMonotonyOffsets(void *data) {
      outputMonotonyOffsets_ = data;
    }
    inline void setOutputOffsets(void *data) {
      outputOffsets_ = data;
    }
    inline void setOutputScalars(void *data) {
      outputScalars_ = data;
    }
    inline void setDeltaApproximate(double data) {
      approxT_.setDelta(data);
    }

  protected:
    bool IgnoreBoundary{false};
    ftm::FTMTreePP contourTree_{};
    dcg::DiscreteGradient dcg_{};
    PersistentSimplexPairs psp_{};
    DiscreteMorseSandwichMPI dms_{};

    // int BackEnd{0};
    BACKEND BackEnd{BACKEND::DISCRETE_MORSE_SANDWICH};
    // progressivity
    ttk::ProgressiveTopology progT_{};
    ttk::ApproximateTopology approxT_{};

    int StartingResolutionLevel{0};
    int StoppingResolutionLevel{-1};
    bool UseTasks{true};
    bool IsResumable{false};
    double TimeLimit{};

    // Approximation
    void *outputScalars_{};
    void *outputOffsets_{};
    void *outputMonotonyOffsets_{};
    double Epsilon;
  };
} // namespace ttk

template <typename scalarType>
int ttk::PersistenceDiagram::computeCTPersistenceDiagram(
  ftm::FTMTreePP &tree,
  const std::vector<
    std::tuple<ttk::SimplexId, ttk::SimplexId, scalarType, bool>> &pairs,
  std::vector<PersistencePair> &diagram) const {

  const ttk::SimplexId numberOfPairs = pairs.size();
  diagram.resize(numberOfPairs);
  for(ttk::SimplexId i = 0; i < numberOfPairs; ++i) {
    const ttk::SimplexId v0 = std::get<0>(pairs[i]);
    const ttk::SimplexId v1 = std::get<1>(pairs[i]);
    const bool type = std::get<3>(pairs[i]);

    if(type == true) {
      diagram[i] = PersistencePair{
        CriticalVertex{
          v0,
          {},
          {},
          {},
          getNodeType(tree.getJoinTree(), ftm::TreeType::Join, v0)},
        CriticalVertex{
          v1,
          {},
          {},
          {},
          getNodeType(tree.getJoinTree(), ftm::TreeType::Join, v1)},
        0, true};
    } else {
      diagram[i] = PersistencePair{
        CriticalVertex{
          v1,
          {},
          {},
          {},
          getNodeType(tree.getSplitTree(), ftm::TreeType::Split, v1)},
        CriticalVertex{
          v0,
          {},
          {},
          {},
          getNodeType(tree.getSplitTree(), ftm::TreeType::Split, v0)},
        2, true};
    }
  }

  diagram.back().isFinite = false; // global extrema pair is infinite

  return 0;
}

template <typename scalarType, typename triangulationType>
void ttk::PersistenceDiagram::augmentPersistenceDiagram(
  std::vector<PersistencePair> &persistencePairs,
  const scalarType *const scalars,
  const triangulationType *triangulation) {

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(std::size_t i = 0; i < persistencePairs.size(); ++i) {
    auto &pair{persistencePairs[i]};
    triangulation->getVertexPoint(pair.birth.id, pair.birth.coords[0],
                                  pair.birth.coords[1], pair.birth.coords[2]);
    pair.birth.sfValue = scalars[pair.birth.id];
    triangulation->getVertexPoint(pair.death.id, pair.death.coords[0],
                                  pair.death.coords[1], pair.death.coords[2]);
    pair.death.sfValue = scalars[pair.death.id];
  }
}

template <typename scalarType, class triangulationType>
int ttk::PersistenceDiagram::execute(std::vector<PersistencePair> &CTDiagram,
                                     const scalarType *inputScalars,
                                     const size_t scalarsMTime,
                                     const SimplexId *inputOffsets,
                                     const triangulationType *triangulation) {

  printMsg(ttk::debug::Separator::L1);

  checkProgressivityRequirement(triangulation);
  checkManifold(triangulation);

  Timer tm{};

  switch(BackEnd) {
    case BACKEND::PERSISTENT_SIMPLEX:
      executePersistentSimplex(CTDiagram, inputOffsets, triangulation);
      break;
    case BACKEND::DISCRETE_MORSE_SANDWICH:
      executeDiscreteMorseSandwich(
        CTDiagram, inputScalars, scalarsMTime, inputOffsets, triangulation);
      break;
    case BACKEND::PROGRESSIVE_TOPOLOGY:
      executeProgressiveTopology(CTDiagram, inputOffsets, triangulation);
      break;
    case BACKEND::APPROXIMATE_TOPOLOGY:
      executeApproximateTopology(CTDiagram, inputScalars, triangulation);
      break;
    case BACKEND::FTM:
      executeFTM(CTDiagram, inputScalars, inputOffsets, triangulation);
      break;
    default:
      printErr("No method was selected");
  }

  this->printMsg("Complete", 1.0, tm.getElapsedTime(), this->threadNumber_);

  if(!isRunningWithMPI()) {
    // augment persistence pairs with meta-data
    augmentPersistenceDiagram(CTDiagram, inputScalars, triangulation);

    // finally sort the diagram
    sortPersistenceDiagram(CTDiagram, inputOffsets);
  } else {
    /* TODO: nothing?*/
  }

  printMsg(ttk::debug::Separator::L1);

  return 0;
}

template <class triangulationType>
int ttk::PersistenceDiagram::executePersistentSimplex(
  std::vector<PersistencePair> &CTDiagram,
  const SimplexId *inputOffsets,
  const triangulationType *triangulation) {

  Timer const tm{};
  const auto dim = triangulation->getDimensionality();

  std::vector<ttk::PersistentSimplexPairs::PersistencePair> pairs{};

  psp_.setDebugLevel(this->debugLevel_);
  psp_.setThreadNumber(this->threadNumber_);
  psp_.computePersistencePairs(pairs, inputOffsets, *triangulation);
  dms_.setInputOffsets(inputOffsets);

  // convert PersistentSimplex pairs (with critical cells id) to PL
  // pairs (with vertices id)

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < pairs.size(); ++i) {
    auto &pair{pairs[i]};
    if(pair.type > 0) {
      pair.birth = dms_.getCellGreaterVertex(
        Cell{pair.type, pair.birth}, *triangulation);
    }
    if(pair.death != -1) {
      pair.death = dms_.getCellGreaterVertex(
        Cell{pair.type + 1, pair.death}, *triangulation);
    }
  }

  CTDiagram.reserve(pairs.size() + 1);

  // find the global maximum
  const auto nVerts = triangulation->getNumberOfVertices();
  const SimplexId globmax = std::distance(
    inputOffsets, std::max_element(inputOffsets, inputOffsets + nVerts));

  // convert pairs to the relevant format
  for(const auto &p : pairs) {
    const auto isFinite = (p.death >= 0);
    const auto death = isFinite ? p.death : globmax;
    if(p.type == 0) {
      const auto deathType = (isFinite && dim > 1)
                               ? CriticalType::Saddle1
                               : CriticalType::Local_maximum;
      CTDiagram.emplace_back(PersistencePair{
        CriticalVertex{p.birth, {}, {}, {}, CriticalType::Local_minimum},
        CriticalVertex{death, {}, {}, {}, deathType}, p.type, isFinite});
    } else if(p.type == 1) {
      const auto birthType
        = (dim == 3) ? CriticalType::Saddle1 : CriticalType::Saddle2;
      const auto deathType = (isFinite && dim == 3)
                               ? CriticalType::Saddle2
                               : CriticalType::Local_maximum;
      CTDiagram.emplace_back(PersistencePair{
        CriticalVertex{p.birth, {}, {}, {}, birthType},
        CriticalVertex{death, {}, {}, {}, deathType}, p.type, isFinite});
    } else if(p.type == 2) {
      CTDiagram.emplace_back(PersistencePair{
        CriticalVertex{p.birth, {}, {}, {}, CriticalType::Saddle2},
        CriticalVertex{death, {}, {}, {}, CriticalType::Local_maximum}, p.type,
        isFinite});
    }
  }

  return 0;
}

template <typename scalarType, class triangulationType>
int ttk::PersistenceDiagram::executeDiscreteMorseSandwich(
  std::vector<PersistencePair> &CTDiagram,
  const scalarType *inputScalars,
  const size_t scalarsMTime,
  const SimplexId *inputOffsets,
  const triangulationType *triangulation) {

  Timer const tm{};
  const auto dim = triangulation->getDimensionality();
  dms_.setUseTasks(UseTasks);

  dms_.buildGradient(inputScalars, scalarsMTime, inputOffsets, *triangulation);
  std::vector<DiscreteMorseSandwichMPI::PersistencePair> dms_pairs{};
  dms_.computePersistencePairs(
    dms_pairs, inputOffsets, *triangulation, this->IgnoreBoundary);
  CTDiagram.resize(dms_pairs.size());

  // transform DiscreteMorseSandwich pairs (critical cells id) to PL
  // pairs (vertices id)
  struct dataRequest {
    ttk::SimplexId gid_;
    ttk::SimplexId lid_;
    char dim_;
    char isBirth_;
  };
  // find the global maximum
  const auto nVerts = triangulation->getNumberOfVertices();
  const SimplexId localMaxId = std::distance(
    inputOffsets, std::max_element(inputOffsets, inputOffsets + nVerts));
  MPI_Datatype MPI_SimplexId = getMPIType(localMaxId);
  struct {
    long int offset;
    int rank;
  } localMax, globalMax;
  localMax.rank = triangulation->getVertexRank(localMaxId);
  localMax.offset = static_cast<long int>(inputOffsets[localMaxId]);
  MPI_Allreduce(
    &localMax, &globalMax, 1, MPI_LONG_INT, MPI_MAXLOC, ttk::MPIcomm_);
  ttk::SimplexId globmax{-1};
  double maxMetaData[4];
  if(globalMax.rank == ttk::MPIrank_) {
    float coords[3];
    globmax = triangulation->getVertexGlobalId(localMaxId);
    triangulation->getVertexPoint(localMaxId, coords[0], coords[1], coords[2]);
    maxMetaData[0] = coords[0];
    maxMetaData[1] = coords[1];
    maxMetaData[2] = coords[2];
    maxMetaData[3] = inputScalars[localMaxId];
  }
  MPI_Bcast(&globmax, 1, MPI_SimplexId, globalMax.rank, ttk::MPIcomm_);
  MPI_Bcast(maxMetaData, 4, MPI_DOUBLE, globalMax.rank, ttk::MPIcomm_);

  std::vector<std::vector<dataRequest>> sendRecvBuffer(ttk::MPIsize_);

  const auto createDataRequestMPIType =
    [this, &MPI_SimplexId](MPI_Datatype &MPI_MessageType) {
      MPI_Datatype types[] = {MPI_SimplexId, MPI_SimplexId, MPI_CHAR, MPI_CHAR};
      int lengths[] = {1, 1, 1, 1};
      const long int mpi_offsets[]
        = {offsetof(dataRequest, gid_), offsetof(dataRequest, lid_),
           offsetof(dataRequest, dim_), offsetof(dataRequest, isBirth_)};
      MPI_Type_create_struct(4, lengths, mpi_offsets, types, &MPI_MessageType);
      MPI_Type_commit(&MPI_MessageType);
    };

  struct dataResponse {
    ttk::SimplexId lid_{-1};
    ttk::SimplexId vertexGid_{-1};
    ttk::SimplexId offset_{-1};
    float coords_[3];
    double sfValue_{0};
    char isBirth_{0};
  };

  const auto createDataResponseMPIType =
    [this, &MPI_SimplexId](MPI_Datatype &MPI_MessageType) {
      MPI_Datatype types[] = {MPI_SimplexId, MPI_SimplexId, MPI_SimplexId,
                              MPI_FLOAT,     MPI_DOUBLE,    MPI_CHAR};
      int lengths[] = {1, 1, 1, 3, 1, 1};
      const long int mpi_offsets[]
        = {offsetof(dataResponse, lid_),     offsetof(dataResponse, vertexGid_),
           offsetof(dataResponse, offset_),  offsetof(dataResponse, coords_),
           offsetof(dataResponse, sfValue_), offsetof(dataResponse, isBirth_)};
      MPI_Type_create_struct(6, lengths, mpi_offsets, types, &MPI_MessageType);
      MPI_Type_commit(&MPI_MessageType);
    };

  const auto fillBirthData
    = [this, &dim](PersistencePair &CTPair,
                   DiscreteMorseSandwichMPI::PersistencePair &p,
                   ttk::SimplexId birthId) {
        CTPair.birth.id = birthId;
        if(p.type == 0) {
          CTPair.birth.type = CriticalType::Local_minimum;
        } else if(p.type == 1) {
          CTPair.birth.type
            = (dim == 3) ? CriticalType::Saddle1 : CriticalType::Saddle2;
        } else if(p.type == 2) {
          CTPair.birth.type = ((p.death >= 0) || dim == 3)
                                ? CriticalType::Saddle2
                                : CriticalType::Local_maximum;
        }
      };

  const auto augmentBirthPersistence
    = [this, &triangulation, &inputOffsets](PersistencePair &CTPair,
                                            ttk::SimplexId lid,
                                            const scalarType *scalars) {
        triangulation->getVertexPoint(lid, CTPair.birth.coords[0],
                                      CTPair.birth.coords[1],
                                      CTPair.birth.coords[2]);
        CTPair.birth.sfValue = scalars[lid];
        CTPair.birth.offset = inputOffsets[lid];
      };

  const auto fillDeathData = [this, &dim](
                               PersistencePair &CTPair,
                               DiscreteMorseSandwichMPI::PersistencePair &p,
                               ttk::SimplexId deathId) {
    const auto isFinite = (p.death >= 0);
    CTPair.death.id = deathId;
    if(p.type == 0) {
      CTPair.death.type = (isFinite && dim > 1) ? CriticalType::Saddle1
                                                : CriticalType::Local_maximum;
    } else if(p.type == 1) {
      CTPair.death.type = (isFinite && dim == 3) ? CriticalType::Saddle2
                                                 : CriticalType::Local_maximum;
    } else if(p.type == 2) {
      CTPair.death.type = CriticalType::Local_maximum;
    }
  };

  const auto augmentDeathPersistence
    = [this, &triangulation, &inputOffsets](PersistencePair &CTPair,
                                            ttk::SimplexId lid,
                                            const scalarType *scalars) {
        triangulation->getVertexPoint(lid, CTPair.death.coords[0],
                                      CTPair.death.coords[1],
                                      CTPair.death.coords[2]);
        CTPair.death.sfValue = scalars[lid];
        CTPair.death.offset = inputOffsets[lid];
      };
  const auto getBirthSimplexType = [this, &dim](const int type) {
    switch(type) {
      case 0:
        return 0;
      case 1:
        return 1;
      default:
        if(dim == 3) {
          return 2;
        } else {
          return 1;
        }
    }
  };

  const auto getDeathSimplexType = [this, &dim](const int type) {
    switch(type) {
      case 0:
        return 1;
      case 1:
        return 2;
      default:
        if(dim == 3) {
          return 3;
        } else {
          return 2;
        }
    }
  };
#ifdef TTK_ENABLE_OPENMP
//#pragma omp declare reduction(merge :std::vector<dataRequest>:
//omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end())) #pragma omp
//parallel for num_threads(threadNumber_) reduction(merge:
//sendRecvBuffer.at(ttk::MPIrank_))
#endif // TTK_ENABLE_OPENMP
  MPI_Barrier(ttk::MPIcomm_);
  for(ttk::SimplexId i = 0; i < dms_pairs.size(); ++i) {
    auto &pair{dms_pairs[i]};
    int simplexType = getBirthSimplexType(pair.type);
    ttk::SimplexId lid
      = triangulation->getSimplexLocalId(pair.birth, simplexType);
    if(lid != -1
       && triangulation->getSimplexRank(lid, simplexType) == ttk::MPIrank_) {
      if(pair.type > 0) {
        lid = dms_.getCellGreaterVertex(Cell{pair.type, lid}, *triangulation);
        pair.birth = triangulation->getVertexGlobalId(lid);
      }
      CTDiagram[i].dim = pair.type;
      CTDiagram[i].isFinite = (pair.death >= 0);
      // Add all the other stuff
      fillBirthData(CTDiagram[i], pair, pair.birth);
      augmentBirthPersistence(CTDiagram[i], lid, inputScalars);
    } else {
      sendRecvBuffer[ttk::MPIrank_].emplace_back(
        dataRequest{pair.birth, i, pair.type, 1});
    }
    if(pair.death == -1) {
      CTDiagram[i].dim = pair.type;
      CTDiagram[i].isFinite = (pair.death >= 0);
      pair.death = globmax;
      fillDeathData(CTDiagram[i], pair, pair.death);
      CTDiagram[i].death.coords[0] = maxMetaData[0];
      CTDiagram[i].death.coords[1] = maxMetaData[1];
      CTDiagram[i].death.coords[2] = maxMetaData[2];
      CTDiagram[i].death.sfValue = maxMetaData[3];
      CTDiagram[i].death.offset = globalMax.offset;
    } else {
      simplexType = getDeathSimplexType(pair.type);
      lid = triangulation->getSimplexLocalId(pair.death, simplexType);
      if(lid != -1
         && triangulation->getSimplexRank(lid, simplexType) == ttk::MPIrank_) {
        CTDiagram[i].dim = pair.type;
        CTDiagram[i].isFinite = (pair.death >= 0);
        lid
          = dms_.getCellGreaterVertex(Cell{pair.type + 1, lid}, *triangulation);
        pair.death = triangulation->getVertexGlobalId(lid);
        fillDeathData(CTDiagram[i], pair, pair.death);
        augmentDeathPersistence(CTDiagram[i], lid, inputScalars);
      } else {
        sendRecvBuffer[ttk::MPIrank_].emplace_back(
          dataRequest{pair.death, i, pair.type, 0});
      }
    }
  }
  // Broadcast all the data to everyone
  // First, broadcast the size of the data to send
  std::vector<int> recvBufferSize(ttk::MPIsize_, 0);
  recvBufferSize[ttk::MPIrank_] = sendRecvBuffer[ttk::MPIrank_].size();
  for(int i = 0; i < ttk::MPIsize_; i++) {
    MPI_Bcast(&recvBufferSize[i], 1, MPI_INTEGER, i, ttk::MPIcomm_);
    if(i != ttk::MPIrank_) {
      sendRecvBuffer[i].resize(recvBufferSize[i]);
    }
  }
  // Then, broadcast the actual data
  MPI_Datatype MPI_requestDataType;
  createDataRequestMPIType(MPI_requestDataType);
  for(int i = 0; i < ttk::MPIsize_; i++) {
    MPI_Bcast(sendRecvBuffer[i].data(), recvBufferSize[i], MPI_requestDataType,
              i, ttk::MPIcomm_);
  }
  std::vector<std::vector<dataResponse>> response(ttk::MPIsize_);

  // For each received element:
  // if it is own by the triangulation, get the data and place to send back
  //#pragma omp parallel for schedule(dynamic, 1)
  for(int i = 0; i < ttk::MPIsize_; i++) {
    if(i != ttk::MPIrank_) {
      for(int j = 0; j < recvBufferSize[i]; j++) {
        auto &element{sendRecvBuffer[i][j]};
        int simplexType;
        if(element.isBirth_) {
          simplexType = getBirthSimplexType(element.dim_);
        } else {
          simplexType = getDeathSimplexType(element.dim_);
        }
        ttk::SimplexId lid
          = triangulation->getSimplexLocalId(element.gid_, simplexType);

        if(lid != -1
           && triangulation->getSimplexRank(lid, simplexType)
                == ttk::MPIrank_) {
          // Add the relevant data
          struct dataResponse res {
            .lid_ = element.lid_, .isBirth_ = element.isBirth_
          };
          ttk::SimplexId vLid = dms_.getCellGreaterVertex(
            Cell{element.dim_ + (1 - element.isBirth_), lid}, *triangulation);
          res.vertexGid_ = triangulation->getVertexGlobalId(vLid);
          res.offset_ = inputOffsets[vLid];
          triangulation->getVertexPoint(
            vLid, res.coords_[0], res.coords_[1], res.coords_[2]);
          res.sfValue_ = inputScalars[vLid];
          // Store to send
          response[i].emplace_back(res);
        }
      }
    }
  }
  // Send back the data
  // First, exchange the size of the data to exchange
  std::vector<dataResponse> responseBuffer;
  std::vector<int> sendBufferSize(ttk::MPIsize_, 0);
  std::vector<int> sendDispls(ttk::MPIsize_, 0);
  std::vector<int> recvDispls(ttk::MPIsize_, 0);
  std::vector<dataResponse> recvBuffer;

  sendBufferSize[0] = response[0].size();
  responseBuffer.insert(
    responseBuffer.end(), response[0].begin(), response[0].end());
  for(int i = 1; i < ttk::MPIsize_; i++) {
    sendBufferSize[i] = response[i].size();
    sendDispls[i] = sendDispls[i - 1] + sendBufferSize[i - 1];
    responseBuffer.insert(
      responseBuffer.end(), response[i].begin(), response[i].end());
  }
  MPI_Alltoall(sendBufferSize.data(), 1, MPI_INTEGER, recvBufferSize.data(), 1,
               MPI_INTEGER, ttk::MPIcomm_);
  for(int i = 1; i < ttk::MPIsize_; i++) {
    recvDispls[i] = recvDispls[i - 1] + recvBufferSize[i - 1];
  }
  recvBuffer.resize(recvDispls.back() + recvBufferSize.back());
  MPI_Datatype MPI_responseDataType;
  createDataResponseMPIType(MPI_responseDataType);

  // Then, exchange to actual data
  MPI_Alltoallv(responseBuffer.data(), sendBufferSize.data(), sendDispls.data(),
                MPI_responseDataType, recvBuffer.data(), recvBufferSize.data(),
                recvDispls.data(), MPI_responseDataType, ttk::MPIcomm_);
  // Receive the data and store it appropriately
  for(const auto &element : recvBuffer) {
    auto &CTPair{CTDiagram[element.lid_]};
    if(element.isBirth_) {
      fillBirthData(CTPair, dms_pairs[element.lid_], element.vertexGid_);
      CTPair.birth.coords[0] = element.coords_[0];
      CTPair.birth.coords[1] = element.coords_[1];
      CTPair.birth.coords[2] = element.coords_[2];
      CTPair.birth.sfValue = element.sfValue_;
      CTPair.birth.offset = element.offset_;
    } else {
      fillDeathData(CTPair, dms_pairs[element.lid_], element.vertexGid_);
      CTPair.death.coords[0] = element.coords_[0];
      CTPair.death.coords[1] = element.coords_[1];
      CTPair.death.coords[2] = element.coords_[2];
      CTPair.death.sfValue = element.sfValue_;
      CTPair.death.offset = element.offset_;
    }
  }
  // Create MPI type for Critical Vertex
  const auto createCriticalVertexMPIType =
    [this, &MPI_SimplexId](MPI_Datatype &MPI_MessageType) {
      MPI_Datatype types[]
        = {MPI_SimplexId, MPI_DOUBLE, MPI_SimplexId, MPI_FLOAT, MPI_INTEGER};
      int lengths[] = {1, 1, 1, 3, 1};
      const long int mpi_offsets[]
        = {offsetof(CriticalVertex, id), offsetof(CriticalVertex, sfValue),
           offsetof(CriticalVertex, offset), offsetof(CriticalVertex, coords),
           offsetof(CriticalVertex, type)};
      MPI_Type_create_struct(5, lengths, mpi_offsets, types, &MPI_MessageType);
      // printMsg("create_struct done");
      MPI_Type_commit(&MPI_MessageType);
    };
  MPI_Datatype MPI_CriticalVertex;
  createCriticalVertexMPIType(MPI_CriticalVertex);
  // Create MPI type for PersistencePair
  const auto createPersistencePairMPIType =
    [this, &MPI_CriticalVertex, &MPI_SimplexId](MPI_Datatype &MPI_MessageType) {
      MPI_Datatype types[]
        = {MPI_CriticalVertex, MPI_CriticalVertex, MPI_SimplexId, MPI_CHAR};
      int lengths[] = {1, 1, 1, 1};
      const long int mpi_offsets[]
        = {offsetof(PersistencePair, birth), offsetof(PersistencePair, death),
           offsetof(PersistencePair, dim), offsetof(PersistencePair, isFinite)};
      MPI_Type_create_struct(4, lengths, mpi_offsets, types, &MPI_MessageType);
      MPI_Type_commit(&MPI_MessageType);
    };
  MPI_Datatype MPI_PersistencePair;
  createPersistencePairMPIType(MPI_PersistencePair);
  // Sort in parallel using the offset of the birth and psort
  std::vector<ttk::SimplexId> vertexDistribution(ttk::MPIsize_);
  ttk::SimplexId localVertexNumber = CTDiagram.size();
  MPI_Allgather(&localVertexNumber, 1, MPI_SimplexId, vertexDistribution.data(),
                1, MPI_SimplexId, ttk::MPIcomm_);
  p_sort::parallel_sort<PersistencePair>(
    CTDiagram, persistenceSort::comp, persistenceSort::oppositeComp,
    vertexDistribution, MPI_PersistencePair, MPI_SimplexId, threadNumber_);
  return 0;
};

template <typename scalarType, class triangulationType>
int ttk::PersistenceDiagram::executeApproximateTopology(
  std::vector<PersistencePair> &CTDiagram,
  const scalarType *inputScalars,
  const triangulationType *triangulation) {

  approxT_.setDebugLevel(debugLevel_);
  approxT_.setThreadNumber(threadNumber_);
  approxT_.setupTriangulation(const_cast<ttk::ImplicitTriangulation *>(
    (const ImplicitTriangulation *)triangulation));
  approxT_.setStartingResolutionLevel(StartingResolutionLevel);
  approxT_.setStoppingResolutionLevel(StoppingResolutionLevel);
  approxT_.setPreallocateMemory(true);
  approxT_.setEpsilon(Epsilon);

  std::vector<ApproximateTopology::PersistencePair> resultDiagram{};

  approxT_.computeApproximatePD(
    resultDiagram, inputScalars, (scalarType *)outputScalars_,
    (SimplexId *)outputOffsets_, (int *)outputMonotonyOffsets_);

  // create the final diagram
  for(const auto &p : resultDiagram) {
    if(p.pairType == 0) {
      CTDiagram.emplace_back(PersistencePair{
        CriticalVertex{p.birth, {}, {}, {}, CriticalType::Local_minimum},
        CriticalVertex{p.death, {}, {}, {}, CriticalType::Saddle1}, p.pairType,
        true});
    } else if(p.pairType == 2) {
      CTDiagram.emplace_back(PersistencePair{
        CriticalVertex{p.birth, {}, {}, {}, CriticalType::Saddle2},
        CriticalVertex{p.death, {}, {}, {}, CriticalType::Local_maximum},
        p.pairType, true});
    } else if(p.pairType == -1) {
      CTDiagram.emplace_back(PersistencePair{
        CriticalVertex{p.birth, {}, {}, {}, CriticalType::Local_minimum},
        CriticalVertex{p.death, {}, {}, {}, CriticalType::Local_maximum},
        p.pairType, false});
    }
  }

  return 0;
}

template <class triangulationType>
int ttk::PersistenceDiagram::executeProgressiveTopology(
  std::vector<PersistencePair> &CTDiagram,
  const SimplexId *inputOffsets,
  const triangulationType *triangulation) {

  progT_.setDebugLevel(debugLevel_);
  progT_.setThreadNumber(threadNumber_);
  progT_.setupTriangulation(const_cast<ttk::ImplicitTriangulation *>(
    (const ImplicitTriangulation *)triangulation));
  progT_.setStartingResolutionLevel(StartingResolutionLevel);
  progT_.setStoppingResolutionLevel(StoppingResolutionLevel);
  progT_.setTimeLimit(TimeLimit);
  progT_.setIsResumable(IsResumable);
  progT_.setPreallocateMemory(true);

  std::vector<ProgressiveTopology::PersistencePair> resultDiagram{};

  progT_.computeProgressivePD(resultDiagram, inputOffsets);

  // create the final diagram
  for(const auto &p : resultDiagram) {
    if(p.pairType == 0) {
      CTDiagram.emplace_back(PersistencePair{
        CriticalVertex{p.birth, {}, {}, {}, CriticalType::Local_minimum},
        CriticalVertex{p.death, {}, {}, {}, CriticalType::Saddle1}, p.pairType,
        true});
    } else if(p.pairType == 2) {
      CTDiagram.emplace_back(PersistencePair{
        CriticalVertex{p.birth, {}, {}, {}, CriticalType::Saddle2},
        CriticalVertex{p.death, {}, {}, {}, CriticalType::Local_maximum},
        p.pairType, true});
    } else if(p.pairType == -1) {
      CTDiagram.emplace_back(PersistencePair{
        CriticalVertex{p.birth, {}, {}, {}, CriticalType::Local_minimum},
        CriticalVertex{p.death, {}, {}, {}, CriticalType::Local_maximum}, 0,
        false});
    }
  }

  return 0;
}

template <typename scalarType, class triangulationType>
int ttk::PersistenceDiagram::executeFTM(
  std::vector<PersistencePair> &CTDiagram,
  const scalarType *inputScalars,
  const SimplexId *inputOffsets,
  const triangulationType *triangulation) {

  contourTree_.setVertexScalars(inputScalars);
  contourTree_.setTreeType(ftm::TreeType::Join_Split);
  contourTree_.setVertexSoSoffsets(inputOffsets);
  contourTree_.setSegmentation(false);
  contourTree_.build<scalarType>(triangulation);

  // get persistence pairs
  std::vector<std::tuple<ttk::SimplexId, ttk::SimplexId, scalarType>> JTPairs;
  std::vector<std::tuple<ttk::SimplexId, ttk::SimplexId, scalarType>> STPairs;
  contourTree_.computePersistencePairs<scalarType>(JTPairs, true);
  contourTree_.computePersistencePairs<scalarType>(STPairs, false);

  // merge pairs
  const auto JTSize = JTPairs.size();
  const auto STSize = STPairs.size();
  std::vector<std::tuple<ttk::SimplexId, ttk::SimplexId, scalarType, bool>>
    CTPairs(JTSize + STSize);
  for(size_t i = 0; i < JTSize; ++i) {
    const auto &x = JTPairs[i];
    CTPairs[i]
      = std::make_tuple(std::get<0>(x), std::get<1>(x), std::get<2>(x), true);
  }
  for(size_t i = 0; i < STSize; ++i) {
    const auto &x = STPairs[i];
    CTPairs[JTSize + i]
      = std::make_tuple(std::get<0>(x), std::get<1>(x), std::get<2>(x), false);
  }

  // remove the last pair which is present two times (global extrema pair)
  if(!CTPairs.empty()) {
    auto cmp =
      [](
        const std::tuple<ttk::SimplexId, ttk::SimplexId, scalarType, bool> &a,
        const std::tuple<ttk::SimplexId, ttk::SimplexId, scalarType, bool> &b) {
        return std::get<2>(a) < std::get<2>(b);
      };

    std::sort(CTPairs.begin(), CTPairs.end(), cmp);
    CTPairs.erase(CTPairs.end() - 1);
  }

  // get persistence diagrams
  computeCTPersistenceDiagram<scalarType>(contourTree_, CTPairs, CTDiagram);

  return 0;
}

template <class triangulationType>
void ttk::PersistenceDiagram::checkProgressivityRequirement(
  const triangulationType *ttkNotUsed(triangulation)) {

  if((BackEnd == BACKEND::PROGRESSIVE_TOPOLOGY
      || BackEnd == BACKEND::APPROXIMATE_TOPOLOGY)
     && !std::is_same<ttk::ImplicitWithPreconditions, triangulationType>::value
     && !std::is_same<ttk::ImplicitNoPreconditions, triangulationType>::value) {

    printWrn("Explicit, Compact or Periodic triangulation detected.");
    printWrn("Defaulting to the FTM backend.");

    BackEnd = BACKEND::FTM;
  }
}

template <class triangulationType>
void ttk::PersistenceDiagram::checkManifold(
  const triangulationType *const triangulation) {

  if(this->BackEnd != BACKEND::DISCRETE_MORSE_SANDWICH) {
    return;
  }

  if(!triangulation->isManifold()) {
    this->printWrn("Non-manifold data-set detected.");
    this->printWrn("Defaulting to the Persistence Simplex backend.");

    this->BackEnd = BACKEND::PERSISTENT_SIMPLEX;
  }
}
