#ifndef DUNE_MULTIDOMAINGRID_GMSHREADERPARSER_HH
#define DUNE_MULTIDOMAINGRID_GMSHREADERPARSER_HH

#include <dune/common/exceptions.hh>

namespace Dune {

namespace mdgrid {

// forward declaration
template<typename HostGrid,typename MDGridTraits>
class MultiDomainGrid;

} // namespace mdgrid

// forward declaration
template<typename GridType>
class GmshReaderParser;

template<class HostGrid, class MDGTraits>
class GmshReaderParser<Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>> : public Dune::GmshReaderParser<HostGrid> {
public:

  using Grid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;

  GmshReaderParser(Dune::GridFactory<Grid>& factory, bool v, bool i)
    : Dune::GmshReaderParser<HostGrid>{factory.hostGridFactory(), v, i}
    , _factory{factory}
  {}

  void read (const std::string& f)
  {
    Dune::GmshReaderParser<HostGrid>::read(f);
    // list of gmsh domain indices for each entity (gmsh index starts at 1!)
    const auto& sub_domain_ids = this->elementIndexMap();
    if (begin(sub_domain_ids) == end(sub_domain_ids))
      return;
    int max_subdomains =
      *std::max_element(begin(sub_domain_ids), end(sub_domain_ids));
    if (max_subdomains < 1)
      DUNE_THROW(IOError,
                 "Gmsh phisical entity indices are 1-index based, but a 0 "
                 "index was read");
    _factory.makeGrid(max_subdomains);
    auto& grid = _factory.grid();

    grid.startSubDomainMarking();
    unsigned int i = 0;
    for (const auto& cell : elements(grid.leafGridView())) {
      assert(sub_domain_ids[i] != 0);
      grid.addToSubDomain(sub_domain_ids[i] - 1, cell), i++;
    }

    grid.preUpdateSubDomains();
    grid.updateSubDomains();
    grid.postUpdateSubDomains();
  }

private:
  Dune::GridFactory<Grid>& _factory;
};

} // namespace Dune

#endif // DUNE_MULTIDOMAINGRID_GMSHREADERPARSER_HH
