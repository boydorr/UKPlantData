
cd("/home/claireh/Documents/UK")
using JuliaDB
using JuliaDBMeta

colnames = ["NBN_Atlas_record_ID", "Occurrence_ID", "Licence",	"Rightsholder",	"Scientific_name", "Taxon_author",	"Name_qualifier", "Common name",	"Species ID (TVK)",	"Taxon Rank",	"Occurrence status",	"Start date",	"Start date day",	"Start date month",	"Start date year",	"End date",	"End date day",	"End date month",	"End date year",	"Locality",	"OSGR",	"Latitude",	"Longitude",	"Coordinate uncertainty (m)",	"Verbatim depth",	"Recorder",	"Determiner",	"Individual count",	"Abundance",	"Abundance scale",	"Organism scope",	"Organism remarks",	"Sex",	"Life stage",	"Occurrence remarks",	"Identification verification status",	"Basis of record",	"Survey key",	"Dataset name",	"Dataset ID",	"Data provider",	"Data provider ID",	"Institution code",	"Kingdom",	"Phylum",	"Class",	"Order",	"Family",	"Genus",	"OSGR 100km",	"OSGR 10km",	"OSGR 2km",	"OSGR 1km",	"Country",	"State/Province",	"Vitality"]
coldict = Dict(:Startdate => String, :Enddate => String)
bry = loadtable("Bryophytes/records1-2020-02-10.csv", colparsers = coldict, type_detect_rows = 10000)
bry = @transform bry {Latitude = convert(Float64, :Latitude)}
