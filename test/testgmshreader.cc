
#include "config.h"
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/multidomaingrid.hh>
#include <iostream>

template<class Grid>
void test_gmsh_reader(std::string fileName, std::size_t expectedDomains) {
  auto mdgrid = Dune::GmshReader<Grid>::read(fileName);

  if (mdgrid->maxAssignedSubDomainIndex() + 1 != static_cast<typename Grid::SubDomainIndex>(expectedDomains))
    DUNE_THROW(Dune::InvalidStateException, "Grid does not contain expected number of sub-domains");

  gridinfo(*mdgrid);
  for (typename Grid::SubDomainIndex i = 0; i != mdgrid->maxAssignedSubDomainIndex()+1; ++i)
    gridinfo(mdgrid->subDomain(i));
}

int main(int argc, char** argv)
{

  try {

    Dune::MPIHelper::instance(argc,argv);

    if (argc != 3) {
      std::cerr << "Usage: " << argv[0] << " <gmhs-file> <expected-domains>" << std::endl;
      return 1;
    }

    std::string fileName = argv[1];
    std::size_t expectedDomains = atoi(argv[2]);

    using HostGrid = Dune::UGGrid<2>;

    test_gmsh_reader<Dune::MultiDomainGrid<HostGrid,Dune::mdgrid::FewSubDomainsTraits<2,8> >>(fileName, expectedDomains);
    test_gmsh_reader<Dune::MultiDomainGrid<HostGrid,Dune::mdgrid::DynamicSubDomainCountTraits<2,1> >>(fileName, expectedDomains);

  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

}
