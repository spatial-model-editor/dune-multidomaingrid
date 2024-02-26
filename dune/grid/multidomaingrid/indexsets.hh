#ifndef DUNE_MULTIDOMAINGRID_INDEXSETS_HH
#define DUNE_MULTIDOMAINGRID_INDEXSETS_HH

#include <map>
#include <unordered_map>
#include <vector>
#include <array>
#include <algorithm>
#include <memory>
#include <type_traits>
#include <tuple>
#include <optional>
#include <utility>

#include <dune/common/hybridutilities.hh>

#include <dune/geometry/typeindex.hh>

#include <dune/grid/common/exceptions.hh>
#include <dune/grid/common/indexidset.hh>

#include <dune/grid/multidomaingrid/utility.hh>
#include <dune/grid/multidomaingrid/subdomaingrid/indexsets.hh>

namespace Dune {

namespace mdgrid {

template<typename HostGrid, typename MDGridTraits>
class MultiDomainGrid;

//! @cond DEV_DOC

//! \internal
namespace detail {

template<template<int> class StructForCodim, typename Codims>
struct _buildMap;

//! \internal template meta program for assembling the per-codim map vector
template<template<int> class StructForCodim, std::size_t... codim>
struct _buildMap<StructForCodim,std::index_sequence<codim...> > {

  typedef std::tuple<StructForCodim<codim>... > type;

};

template<template<int> class StructForCodim, int dimension>
struct buildMap {

  typedef typename _buildMap<
    StructForCodim,
    decltype(std::make_index_sequence<dimension+1>())
    >::type type;

};

}

//! @endcond

/**
 * This is the heart of this module. The index sets for each level of the grid view
 * and for the leaf grid view store a container of the type ContainerMap which is
 * roughly the following:
 *
 * @code{c++}
 * template<int codim>
 * struct Container {
 *  struct MapEntry {
 *    using SubDomainSet = typename MDGridTraits::template Codim<cc>::SubDomainSet;
 *    SubDomainSet domains;
 *    IndexType index;
 *  };
 *
 *  // [geometry-type, host-index] -> map-entry
 *  using IndexMap = std::vector<std::vector<MapEntry>>;
 *
 *  // [geometry-type, sub-domain] -> entity-count
 *  using SizeMap = std::vector<typename Grid::MDGridTraits::template Codim<codim>::SizeContainer>;
 *
 *  // [sub-domain] -> entity-count
 *  using CodimSizeMap = typename Grid::MDGridTraits::template Codim<codim>::SizeContainer;
 *
 *  // [map-entry.index] -> [sub-domain-offset-0-index, ..., sub-domain-offset-i-index]
 *  using MultiIndexMap = std::vector<typename Grid::MDGridTraits::template Codim<codim>::MultiIndexContainer>;

 *  IndexMap indexMap;
 *  SizeMap sizeMap;
 *  CodimSizeMap codimSizeMap;
 *  MultiIndexMap multiIndexMap;
 * };
 *
 * using ContainerMap = std::tuple<Container<0>, ..., Container<dim>>;
 * @endcode
 *
 * Most of the fancy stuff is just to avoid instantiating unsupported codimensions by the host grid.
 * But in general, the index of the sub-domain entity is stored in `multiIndexMap` and is
 * accessed through the `indexMap`. If there is only one sub-domain in the entity, the sub-domain
 * index is directly stored in `indexMap` to avoid the indirection of going to `multiIndexMap`.
 * This also makes the multi-domain entities more compact.
 *
 * To be performant in the overlapping case, the current implementation requires that
 * `MultiIndexContainer` stores its contents within its class itself (i.e. no indirection to the heap).
 * Doing otherwise implies very scatter memory and a very irregular memory access pattern. Either a
 * fixed integer or an array do this. An idea to make this performant in a more general case is to use
 * a polymorphic allocator that creates the multi-indices contiguously in memory.
 */

/**
 * Index set for the MultiDomainGrid.
 */
template<typename GridImp, typename HostGridViewType>
class IndexSetWrapper :
    public Dune::IndexSet<GridImp,IndexSetWrapper<GridImp,HostGridViewType>,
                          typename HostGridViewType::IndexSet::IndexType,
                          typename HostGridViewType::IndexSet::Types>
{

  using Grid = std::remove_const_t<GridImp>;

  template<typename, typename>
  friend class IndexSetWrapper;

  template<typename, typename>
  friend class MultiDomainGrid;

  template<typename, typename>
  friend class subdomain::IndexSetWrapper;

  template<typename, typename, typename, typename>
  friend class SubDomainInterface;

  template<typename>
  friend class SubDomainToSubDomainController;

  template<typename>
  friend class AllInterfacesController;

  typedef IndexSetWrapper<GridImp,HostGridViewType> ThisType;

  using HostGrid = typename Grid::HostGrid;
  typedef HostGridViewType HostGridView;
  typedef typename HostGridView::IndexSet HostIndexSet;
  using ctype = typename Grid::ctype;
  using CellReferenceElement = Dune::ReferenceElement<typename HostGrid::template Codim<0>::Entity::Geometry>;

  typedef Dune::IndexSet<
    GridImp,
    IndexSetWrapper<
      GridImp,
      HostGridViewType
      >,
    typename HostGridViewType::IndexSet::IndexType,
    typename HostGridViewType::IndexSet::Types
    > BaseT;

public:

  typedef typename BaseT::Types Types;

  using MDGridTraits = typename Grid::MDGridTraits;
  typedef typename MDGridTraits::template Codim<0>::SubDomainSet SubDomainSet;
  typedef typename MDGridTraits::SubDomainIndex SubDomainIndex;


  typedef typename HostIndexSet::IndexType IndexType;
  static const int dimension = Grid::dimension;
  static const std::size_t maxSubDomains = SubDomainSet::maxSize;

private:

  typedef typename HostGridView::template Codim<0>::Iterator HostEntityIterator;
  typedef typename HostGridView::template Codim<0>::Entity HostEntity;
  using Codim0Entity = typename Grid::Traits::template Codim<0>::Entity;

  template<int cc>
  struct MapEntry {

    typedef typename MDGridTraits::template Codim<cc>::SubDomainSet SubDomainSet;

    SubDomainSet domains;
    IndexType index;
  };

  //! Placeholder struct for non-supported codimensions in data structures.
  struct NotSupported {};

  template<int c>
  struct Containers {
    inline static const int codim = c;

    inline static const bool supported = Grid::MDGridTraits::template Codim<codim>::supported;

    static_assert((codim > 0 || supported), "index mapping of codimension 0 must be supported!");

    using IndexMap = std::conditional_t<
      supported,
      std::vector<std::vector<MapEntry<codim> > >,
      NotSupported
      >;

    using SizeMap = std::conditional_t<
      supported,
      std::vector<typename Grid::MDGridTraits::template Codim<codim>::SizeContainer>,
      NotSupported
      >;

    using CodimSizeMap = std::conditional_t<
      supported,
      typename Grid::MDGridTraits::template Codim<codim>::SizeContainer,
      NotSupported
      >;

    using MultiIndexMap = std::conditional_t<
      supported,
      std::vector<typename Grid::MDGridTraits::template Codim<codim>::MultiIndexContainer>,
      NotSupported
      >;

    IndexMap indexMap;
    SizeMap sizeMap;
    CodimSizeMap codimSizeMap;
    MultiIndexMap multiIndexMap;

    // containers should not be assignable...
    Containers& operator=(const Containers&) = delete;

    // ...but must be movable
    Containers& operator=(Containers&&) = default;

    // make sure the compiler generates all necessary constructors...
    Containers(const Containers&) = default;
    Containers(Containers&&) = default;
    Containers() = default;

  };

  typedef typename detail::buildMap<Containers,dimension>::type ContainerMap;

  typedef std::vector<std::shared_ptr<IndexSetWrapper<GridImp, typename HostGridView::Grid::LevelGridView> > > LevelIndexSets;

public:

  //! Returns the index of the entity with codimension codim.
  template<int codim>
  IndexType index(const typename Grid::Traits::template Codim<codim>::Entity& e) const {
    return _hostGridView.indexSet().index(_grid.hostEntity(e));
  }

  //! Returns the index of the entity.
  template<typename Entity>
  IndexType index(const Entity& e) const {
    return _hostGridView.indexSet().index(_grid.hostEntity(e));
  }

  //! Returns the subdindex of the i-th subentity of e with codimension codim.
  template<int codim>
  IndexType subIndex(const typename Grid::Traits::template Codim<codim>::Entity& e, int i, unsigned int cd) const {
    return _hostGridView.indexSet().template subIndex<codim>(_grid.hostEntity(e),i,cd);
  }

  //! Returns the subdindex of the i-th subentity of e with codimension codim.
  template<typename Entity>
  IndexType subIndex(const Entity& e, int i, unsigned int cd) const {
    return _hostGridView.indexSet().subIndex(_grid.hostEntity(e),i,cd);
  }

    //! Returns a list of all geometry types with codimension codim contained in the grid.
  Types types(int codim) const {
    return _hostGridView.indexSet().types(codim);
  }

  //! Returns the number of entities with GeometryType type in the grid.
  IndexType size(GeometryType type) const {
    return _hostGridView.indexSet().size(type);
  }

  //! Returns the number of entities with codimension codim in the grid.
  IndexType size(int codim) const {
    return _hostGridView.indexSet().size(codim);
  }

  //! Returns true if the entity is contained in the grid.
  template<typename EntityType>
  bool contains(const EntityType& e) const {
    return _hostGridView.indexSet().contains(_grid.hostEntity(e));
  }

  //! Returns a constant reference to the SubDomainSet of the given entity.
  template<typename EntityType>
  const typename MapEntry<EntityType::codimension>::SubDomainSet& subDomains(const EntityType& e) const {
    return subDomainsForHostEntity(_grid.hostEntity(e));
  }

  //! Returns a constant reference to the SubDomainSet of the given entity with codimension cc.
  //! \tparam cc the codimension of the entity.
  template<int cc>
  const typename MapEntry<cc>::SubDomainSet& subDomains(const typename Grid::Traits::template Codim<cc>::Entity& e) const {
    return subDomainsForHostEntity<cc>(_grid.hostEntity(e));
  }

  //! Returns the index of the entity in a specific subdomain.
  template<class EntityType>
  IndexType index(SubDomainIndex subDomain, const EntityType& e) const {
    return index<EntityType::codimension>(subDomain,e);
  }

  //! Returns the index of the entity with codimension cc in a specific subdomain.
  //! \tparam the codimension of the entity.
  template<int cc>
  IndexType index(SubDomainIndex subDomain, const typename Grid::Traits::template Codim<cc>::Entity& e) const {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    const MapEntry<cc>& me = indexMap<cc>()[LocalGeometryTypeIndex::index(gt)][hostIndex];
    assert(me.domains.contains(subDomain));
    if (me.domains.simple()) {
      return me.index;
    } else {
      return multiIndexMap<cc>()[me.index][me.domains.domainOffset(subDomain)];
    }
  }

private:

  //! Returns a mutable reference to the SubDomainSet of the given entity.
  template<typename EntityType>
  typename MapEntry<EntityType::codimension>::SubDomainSet& subDomains(const EntityType& e) {
    return subDomainsForHostEntity(_grid.hostEntity(e));
  }

  //! Returns a mutable reference to the SubDomainSet of the given entity with codimension cc.
  //! \tparam cc the codimension of the entity.
  template<int cc>
  typename MapEntry<cc>::SubDomainSet& subDomains(const typename Grid::Traits::template Codim<cc>::Entity& e) {
    return subDomainsForHostEntity<cc>(_grid.hostEntity(e));
  }

  //! Returns a mutable reference to the SubDomainSet of the given host entity.
  template<typename HostEntity>
  typename MapEntry<HostEntity::codimension>::SubDomainSet& subDomainsForHostEntity(const HostEntity& e) {
    return subDomainsForHostEntity<HostEntity::codimension>(e);
  }

  //! Returns a mutable reference to the SubDomainSet of the given entity with codimension cc.
  //! \tparam cc the codimension of the entity.
  template<int cc>
  typename MapEntry<cc>::SubDomainSet& subDomainsForHostEntity(const typename Grid::HostGrid::Traits::template Codim<cc>::Entity& he) {
    return indexMap<cc>()[LocalGeometryTypeIndex::index(he.type())][_hostGridView.indexSet().index(he)].domains;
  }

  //! Returns a constant reference to the SubDomainSet of the given host entity.
  template<typename HostEntity>
  const typename MapEntry<HostEntity::codimension>::SubDomainSet& subDomainsForHostEntity(const HostEntity& e) const {
    return subDomainsForHostEntity<HostEntity::codimension>(e);
  }

  //! Returns a constant reference to the SubDomainSet of the given entity with codimension cc.
  //! \tparam cc the codimension of the entity.
  template<int cc>
  const typename MapEntry<cc>::SubDomainSet& subDomainsForHostEntity(const typename Grid::HostGrid::Traits::template Codim<cc>::Entity& he) const {
    return indexMap<cc>()[LocalGeometryTypeIndex::index(he.type())][_hostGridView.indexSet().index(he)].domains;
  }

  template<int cc>
  IndexType indexForSubDomain(SubDomainIndex subDomain, const typename Grid::HostGrid::Traits::template Codim<cc>::Entity& he) const {
    const GeometryType gt = he.type();
    const IndexType hostIndex = _hostGridView.indexSet().index(he);
    const MapEntry<cc>& me = indexMap<cc>()[LocalGeometryTypeIndex::index(gt)][hostIndex];
    assert(me.domains.contains(subDomain));
    if (me.domains.simple()) {
      return me.index;
    } else {
      return multiIndexMap<cc>()[me.index][me.domains.domainOffset(subDomain)];
    }
  }

  template<typename HostEntity>
  IndexType subIndexForSubDomain(SubDomainIndex subDomain, const HostEntity& he, int i, int codim) const {
    std::optional<IndexType> index;
    auto gt = referenceElement(he.geometry()).type(i,codim - he.codimension);
    auto hostIndex = _hostGridView.indexSet().subIndex(he,i,codim);
    Dune::Hybrid::switchCases(std::make_index_sequence<dimension+1>{}, codim, [&](auto icodim){
      if constexpr (MDGridTraits::template Codim<icodim>::supported) {
        const MapEntry<icodim>& me = this->indexMap<icodim>()[LocalGeometryTypeIndex::index(gt)][hostIndex];
        assert(me.domains.contains(subDomain));
        if (me.domains.simple()) {
          index = me.index;
        } else {
          index = this->multiIndexMap<icodim>()[me.index][me.domains.domainOffset(subDomain)];
        }
      } else
        DUNE_THROW(GridError,"invalid codimension specified");
    });
    assert(index.has_value());
    return *index;
  }

  Types typesForSubDomain(SubDomainIndex subDomain, int codim) const {
    return types(codim);
  }

  IndexType sizeForSubDomain(SubDomainIndex subDomain, GeometryType type) const {
    std::optional<IndexType> index;
    Dune::Hybrid::switchCases(Dune::range(index_constant<dimension+1>{}), dimension-type.dim(), [&](auto icodim){
      if constexpr (MDGridTraits::template Codim<icodim>::supported)
        index = this->sizeMap<icodim>()[LocalGeometryTypeIndex::index(type)][subDomain];
    });
    return index.value_or(0);
  }

  IndexType sizeForSubDomain(SubDomainIndex subDomain, int codim) const {
    std::optional<IndexType> index;
    Dune::Hybrid::switchCases(Dune::range(index_constant<dimension+1>{}), codim, [&](auto icodim){
      if constexpr (MDGridTraits::template Codim<icodim>::supported)
        index = this->codimSizes<icodim>()[subDomain];
    });
    return index.value_or(0);
  }

  template<typename EntityType>
  bool containsForSubDomain(SubDomainIndex subDomain, const EntityType& he) const {
    const GeometryType gt = he.type();
    const IndexType hostIndex = _hostGridView.indexSet().index(he);
    const MapEntry<EntityType::codimension>& me = indexMap<EntityType::codimension>()[LocalGeometryTypeIndex::index(gt)][hostIndex];
    return me.domains.contains(subDomain);
  }

public:

  template<typename SubDomainEntity>
  IndexType subIndex(SubDomainIndex subDomain, const SubDomainEntity& e, int i, int codim) const {
    return subIndexForSubDomain(subDomain,_grid.hostEntity(e),i,codim);
  }

  Types types(SubDomainIndex subDomain, int codim) const {
    return types(codim);
  }

  IndexType size(SubDomainIndex subDomain, GeometryType type) const {
    return sizeForSubDomain(subDomain,type);
  }

  IndexType size(SubDomainIndex subDomain, int codim) const {
    return sizeForSubDomain(subDomain,codim);
  }

  //! Returns true if the entity is contained in a specific subdomain.
  template<typename EntityType>
  bool contains(SubDomainIndex subDomain, const EntityType& e) const {
    const GeometryType gt = e.type();
    const IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    const MapEntry<EntityType::codimension>& me = indexMap<EntityType::codimension>()[LocalGeometryTypeIndex::index(gt)][hostIndex];
    return me.domains.contains(subDomain);
  }

private:

  const GridImp& _grid;
  HostGridView _hostGridView;
  ContainerMap _containers;

  void swap(ThisType& rhs) {
    assert(&_grid == &rhs._grid);
    std::swap(_containers,rhs._containers);
  }

  void addToSubDomain(SubDomainIndex subDomain, const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap<0>()[LocalGeometryTypeIndex::index(gt)][hostIndex].domains.add(subDomain);
  }

  void removeFromSubDomain(SubDomainIndex subDomain, const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap<0>()[LocalGeometryTypeIndex::index(gt)][hostIndex].domains.remove(subDomain);
  }

  void removeFromAllSubDomains(const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap<0>()[LocalGeometryTypeIndex::index(gt)][hostIndex].domains.clear();
  }

  void assignToSubDomain(SubDomainIndex subDomain, const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap<0>()[LocalGeometryTypeIndex::index(gt)][hostIndex].domains.set(subDomain);
  }

  void addToSubDomains(const typename MDGridTraits::template Codim<0>::SubDomainSet& subDomains, const Codim0Entity& e) {
    GeometryType gt = e.type();
    IndexType hostIndex = _hostGridView.indexSet().index(_grid.hostEntity(e));
    indexMap<0>()[LocalGeometryTypeIndex::index(gt)][hostIndex].domains.addAll(subDomains);
  }

public:

  // just make these public for std::make_shared() access

  IndexSetWrapper(const GridImp& grid, HostGridView hostGridView) :
    _grid(grid),
    _hostGridView(hostGridView)
  {}

private:

  explicit IndexSetWrapper(const ThisType& rhs) :
    _grid(rhs._grid),
    _hostGridView(rhs._hostGridView),
    _containers(rhs._containers)
    {}


  //! Returns the index map for the given codimension.
  //! \tparam cc The requested codimension.
  template<int cc>
  typename Containers<cc>::IndexMap& indexMap() {
    return std::get<cc>(_containers).indexMap;
  }

  template<int cc>
  typename Containers<cc>::SizeMap& sizeMap() {
    return std::get<cc>(_containers).sizeMap;
  }

  template<int cc>
  typename Containers<cc>::CodimSizeMap& codimSizes() {
    return std::get<cc>(_containers).codimSizeMap;
  }

  template<int cc>
  typename Containers<cc>::MultiIndexMap& multiIndexMap() {
    return std::get<cc>(_containers).multiIndexMap;
  }

  template<int cc>
  const typename Containers<cc>::IndexMap& indexMap() const {
    return std::get<cc>(_containers).indexMap;
  }

  template<int cc>
  const typename Containers<cc>::SizeMap& sizeMap() const {
    return std::get<cc>(_containers).sizeMap;
  }

  template<int cc>
  const typename Containers<cc>::CodimSizeMap& codimSizes() const {
    return std::get<cc>(_containers).codimSizeMap;
  }

  template<int cc>
  const typename Containers<cc>::MultiIndexMap& multiIndexMap() const {
    return std::get<cc>(_containers).multiIndexMap;
  }

  template<typename Functor>
  void forEachContainer(Functor func) const {
    // static loop over container tuple
    Hybrid::forEach(Dune::range(std::tuple_size<ContainerMap>{}),
      [&](auto codim){
        if constexpr (MDGridTraits::template Codim<codim>::supported)
          func(std::get<codim>(_containers));
      });
  }

  template<typename Functor>
  void forEachContainer(Functor func) {
    // static loop over container tuple
    Hybrid::forEach(Dune::range(std::tuple_size<ContainerMap>{}),
      [&](auto codim){
        if constexpr (MDGridTraits::template Codim<codim>::supported)
          func(std::get<codim>(_containers));
      });
  }

  void reset(bool full) {
    const HostIndexSet& his = _hostGridView.indexSet();
    if (full) {
      ContainerMap cm = ContainerMap();
      std::swap(_containers,cm);
    }
    forEachContainer([&](auto& c){
      // setup per-codim sizemap
      constexpr int codim = std::decay_t<decltype(c)>::codim;
      _grid._traits.template setupSizeContainer<codim>(c.codimSizeMap);
      c.indexMap.resize(LocalGeometryTypeIndex::size(dimension-codim));
      c.sizeMap.resize(LocalGeometryTypeIndex::size(dimension-codim));
      for (auto gt : his.types(codim)) {
        const auto gt_index = LocalGeometryTypeIndex::index(gt);
        if (full) {
          // resize index map
          c.indexMap[gt_index].resize(his.size(gt));
        } else if (codim > 0) {
          // clear out marked state for codim > 0 (we cannot keep the old
          // state for subentities, as doing so will leave stale entries if
          // elements are removed from a subdomain
          for (auto& mapEntry : c.indexMap[gt_index])
            mapEntry.domains.clear();
        }
        // setup / reset SizeMap counter
        auto& size_map = c.sizeMap[gt_index];
        _grid._traits.template setupSizeContainer<codim>(size_map);
        std::fill(size_map.begin(),size_map.end(),0);
      }
      // clear MultiIndexMap
      c.multiIndexMap.clear();
    });
  }

  struct updatePerCodimSizes {

    template<int codim>
    void operator()(Containers<codim>& c) const {
      // reset size for this codim to zero
      std::fill(c.codimSizeMap.begin(),c.codimSizeMap.end(),0);
      // collect per-geometrytype sizes into codim size structure
      std::for_each(c.sizeMap.begin(),
                    c.sizeMap.end(),
                    util::collect_elementwise<std::plus<IndexType> >(c.codimSizeMap));
    }

  };

  void update(LevelIndexSets& levelIndexSets, bool full) {
    const HostIndexSet& his = _hostGridView.indexSet();
    // reset(full);
    for (auto& levelIndexSet : levelIndexSets) {
      levelIndexSet->reset(full);
    }

    this->communicateSubDomainSelection();

    auto& im = indexMap<0>();
    auto& sm = sizeMap<0>();
    for (const auto& he : elements(_hostGridView)) {
      // get map entry for host entity
      auto geo = he.geometry();
      auto hgt = geo.type();
      const auto hgt_index = LocalGeometryTypeIndex::index(hgt);
      IndexType hostIndex = his.index(he);
      MapEntry<0>& me = im[hgt_index][hostIndex];

      if (_grid.supportLevelIndexSets()) {
        // include subdomain in entity lvl
        IndexType hostLvlIndex = levelIndexSets[he.level()]->_hostGridView.indexSet().index(he);
        levelIndexSets[he.level()]->template indexMap<0>()[hgt_index][hostLvlIndex].domains.addAll(me.domains);
        // propagate subdomain to all ancestors
        markAncestors(levelIndexSets,he,me.domains);
      }
      updateMapEntry(me,sm[hgt_index],multiIndexMap<0>());
      forEachContainer(markSubIndices(he,me.domains,his,geo));
    }

    propagateBorderEntitySubDomains();

    forEachContainer(updateSubIndices(*this));
    forEachContainer(updatePerCodimSizes());
    for(auto& levelIndexSet : levelIndexSets) {
      levelIndexSet->updateLevelIndexSet();
    }
  }


  void updateLevelIndexSet() {
    const HostIndexSet& his = _hostGridView.indexSet();
    typename Containers<0>::IndexMap& im = indexMap<0>();
    typename Containers<0>::SizeMap& sm = sizeMap<0>();

    communicateSubDomainSelection();

    for (const auto& he : elements(_hostGridView)) {
      auto geo = he.geometry();
      const GeometryType hgt = geo.type();
      const auto hgt_index = LocalGeometryTypeIndex::index(hgt);
      IndexType hostIndex = his.index(he);
      MapEntry<0>& me = im[hgt_index][hostIndex];
      updateMapEntry(me,sm[hgt_index],multiIndexMap<0>());
      forEachContainer(markSubIndices(he,me.domains,his,geo));
    }

    propagateBorderEntitySubDomains();

    forEachContainer(updateSubIndices(*this));
    forEachContainer(updatePerCodimSizes());
  }

  template<int codim, typename SizeContainer, typename MultiIndexContainer, typename Alloc>
  void updateMapEntry(MapEntry<codim>& me, SizeContainer& sizes, std::vector<MultiIndexContainer, Alloc>& multiIndexMap) {
    switch (me.domains.state()) {
    case MapEntry<codim>::SubDomainSet::emptySet:
      break;
    case MapEntry<codim>::SubDomainSet::simpleSet:
      me.index = sizes[*me.domains.begin()]++;
      break;
    case MapEntry<codim>::SubDomainSet::multipleSet:
      me.index = multiIndexMap.size();
      MultiIndexContainer mic;
      for (const auto& subdomain : me.domains)
        mic[me.domains.domainOffset(subdomain)] = sizes[subdomain]++;
      multiIndexMap.push_back(std::move(mic));
    }
  }

  template<typename SubDomainSet>
  void markAncestors(LevelIndexSets& levelIndexSets, HostEntity he, const SubDomainSet& domains) {
    while (he.level() > 0) {
      he = he.father();
      SubDomainSet& fatherDomains =
        levelIndexSets[he.level()]->template indexMap<0>()[LocalGeometryTypeIndex::index(he.type())][levelIndexSets[he.level()]->_hostGridView.indexSet().index(he)].domains;
      if (fatherDomains.containsAll(domains))
        break;
      fatherDomains.addAll(domains);
    }
  }

  struct markSubIndices {

    template<int codim>
    void operator()(Containers<codim>& c) const {
      if (codim == 0)
        return;
      const int size = _refEl.size(codim);
      for (int i = 0; i < size; ++i) {
        IndexType hostIndex = _his.subIndex(_he,i,codim);
        GeometryType gt = _refEl.type(i,codim);
        c.indexMap[LocalGeometryTypeIndex::index(gt)][hostIndex].domains.addAll(_domains);
      }
    }

    typedef typename MapEntry<0>::SubDomainSet& DomainSet;

    const HostEntity& _he;
    DomainSet& _domains;
    const HostIndexSet& _his;
    CellReferenceElement _refEl;

    markSubIndices(const HostEntity& he, DomainSet& domains, const HostIndexSet& his, const typename HostEntity::Geometry& geo) :
      _he(he),
      _domains(domains),
      _his(his),
      _refEl(referenceElement(geo))
    {}

  };

  struct updateSubIndices {

    template<int codim>
    void operator()(Containers<codim>& c) const {
      if (codim == 0)
        return;
      for (std::size_t gt_index = 0,
             gt_end = c.indexMap.size();
           gt_index != gt_end;
           ++gt_index)
        for (auto& im_entry : c.indexMap[gt_index])
          _indexSet.updateMapEntry(im_entry,c.sizeMap[gt_index],c.multiIndexMap);
    }

    ThisType& _indexSet;

    updateSubIndices(ThisType& indexSet) :
      _indexSet(indexSet)
    {}
  };

  bool supportsCodim(int codim) const
  {
    auto codim_seq = std::make_index_sequence<dimension + 1>{};
    return Dune::Hybrid::switchCases(codim_seq, codim, [&](auto icodim) {
      return MDGridTraits::template Codim<icodim>::supported && Dune::Capabilities::hasEntity<HostGrid, icodim>::v;
    }, []{
      return false;
    });
  }

  template<typename Impl>
  struct SubDomainSetDataHandleBase
    : public Dune::CommDataHandleIF<Impl,
                                    typename MapEntry<0>::SubDomainSet::DataHandle::DataType
                                    >
  {
    typedef typename MapEntry<0>::SubDomainSet SubDomainSet;
    typedef typename SubDomainSet::DataHandle DataHandle;

    bool fixedSize(int dim, int codim) const
    {
      return DataHandle::fixedSize(dim,codim);
    }

    template<typename Entity>
    std::size_t size(const Entity& e) const
    {
      return MapEntry<Entity::codimension>::SubDomainSet::DataHandle::size(_indexSet.subDomainsForHostEntity(e));
    }

    template<typename MessageBufferImp, typename Entity>
    void gather(MessageBufferImp& buf, const Entity& e) const
    {
      MapEntry<Entity::codimension>::SubDomainSet::DataHandle::gather(buf,_indexSet.subDomainsForHostEntity(e));
    }

    template<typename MessageBufferImp, typename Entity>
    void scatter(MessageBufferImp& buf, const Entity& e, std::size_t n)
    {
      MapEntry<Entity::codimension>::SubDomainSet::DataHandle::scatter(buf,_indexSet.subDomainsForHostEntity(e),n);
    }

    SubDomainSetDataHandleBase(ThisType& indexSet)
      : _indexSet(indexSet)
    {}

    ThisType& _indexSet;

  };

  struct SelectionDataHandle
    : public SubDomainSetDataHandleBase<SelectionDataHandle>
  {

    bool contains(int dim, int codim) const
    {
      return codim == 0;
    }

    SelectionDataHandle(ThisType& indexSet)
      : SubDomainSetDataHandleBase<SelectionDataHandle>(indexSet)
    {}

  };

  struct BorderPropagationDataHandle
    : public SubDomainSetDataHandleBase<BorderPropagationDataHandle>
  {

    bool contains(int dim, int codim) const
    {
      return codim > 0 && this->_indexSet.supportsCodim(codim);
    }

    BorderPropagationDataHandle(ThisType& indexSet)
      : SubDomainSetDataHandleBase<BorderPropagationDataHandle>(indexSet)
    {}

  };


  void communicateSubDomainSelection()
  {
    SelectionDataHandle dh(*this);
    _hostGridView.template communicate<SelectionDataHandle>(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  }

  void propagateBorderEntitySubDomains()
  {
    BorderPropagationDataHandle dh(*this);
    _hostGridView.communicate(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  }

};

} // namespace mdgrid

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_INDEXSETS_HH
