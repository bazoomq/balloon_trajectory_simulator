import pandas as pd


def create_kml(coordinates, output_file):
    kml_template = """<?xml version="1.0" encoding="UTF-8"?>
    <kml xmlns="http://www.opengis.net/kml/2.2">
      <Document>
        <name>GPS Data</name>

        <!-- Define a style for red points -->
        <Style id="redPoint">
          <IconStyle>
            <color>ff0000ff</color> <!-- Red color (in ABGR format) -->
            <scale>1.0</scale>
            <Icon>
              <href>http://maps.google.com/mapfiles/ms/icons/red-dot.png</href>
            </Icon>
          </IconStyle>
        </Style>
        """

    for coordinate in coordinates:
        lat, lon, alt = coordinate
        placemark_template = f"""
        <Placemark>
          <styleUrl>#redPoint</styleUrl>
          <Point><altitudeMode>absolute</altitudeMode><extrude>1</extrude>

            <coordinates>{lon},{lat},{alt}</coordinates>
          </Point>
          <altitudeMode>absolute</altitudeMode>
        </Placemark>
        """
        kml_template += placemark_template

    kml_template += """
      </Document>
    </kml>
    """

    with open(output_file, "w") as f:
        f.write(kml_template)
    return kml_template

coordinates = [
    # Add more coordinates here
(40.23792239194559,44.50879090807187,1112.4751203515598),
]

output_file = "output.kml"
create_kml(coordinates, output_file)
