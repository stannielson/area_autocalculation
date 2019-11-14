"""
-------------------------------------------------------------------------------
Title:          autocalculateacres

Author:         Stanton K. Nielson
                GIS Specialist
                Bureau of Land Management
                Wyoming High Desert District
                snielson@blm.gov
                
Date:           October 24, 2019

Version:        1.01

Description:    This module provides autocalculation object capabilities for
                calculating acreage in a feature class.
-------------------------------------------------------------------------------
"""


import arcpy
import re


class AcreageCalculator(object):


    """
    AcreageCalculator
    ---------------------------------------------------------------------------
    Performs acreage calculation of features based on the extent of a feature
    class overall or locally, based on extents of individual of individual
    features.  Calculation relies on projection of feature geometry to an
    Albers Equal Area Conic projection (with extent-derived standard parallels
    and central meridian) with an automatically selected transformation (where
    necessary).
    ---------------------------------------------------------------------------
    Parameters:
    - feature_class:    an input feature class
    - gcs_wkid (int):   the well-known identifier for the desired geographic
                        coordinate system to produce the best accuracy and
                        precision (default is 104145, the North American Datum
                        of 1983 with 2011 HARN adjustment)
    - local (boolean):  specification of whether acreage for features should
                        derive from the extents of individual features or the
                        overall extent of the feature class (default is False).
                        
    WARNING! Calculation is not meant for extents that exceed 30 degrees of
    latitude.  Use of local extents will also produce longer processing time
    for larger datasets.
    ---------------------------------------------------------------------------
    Usage: AcreageCalculator(feature_class, {gcs_wkid=104145}, {local=False})
    """


    def __init__(self, feature_class, gcs_wkid=104145, local=False):

        self.features = feature_class
        self.gcs = arcpy.SpatialReference(gcs_wkid)
        self.local = local
        self.main_transformation = self.__get_main_transformation()
        self.overall_extent = self.__get_main_extent()
        self.main_spatial_ref = self.__build_spatial_ref(self.overall_extent)
        self.CalculateAcreage()


    def CalculateAcreage(self):

        """
        Instance-accessible method to calculate acreage.
        """
        self.__calculate_acreage()
        return


    def __build_spatial_ref(self, extent):

        """
        Builds spatial reference using North America Albers Equal Area Conic
        (WKID 102008). Method modifies the spatial reference WKT to match
        specified GCS/datum and extent-derived standard parallels and central
        meridian.  Modified WKT is then instantiated as a SpatialReference
        object.
        """
        src = arcpy.SpatialReference(102008)
        src_sec = self.__modify_secant_parameters(src, extent)
        src_ref = src_sec.exportToString()
        revision = list()
        spheroid = self.__get_spheroid_parameters()
        revision.append((
            self.__find('\'GRS_1980\',[-+]?\d*[.]\d*,[-+]?\d*[.]\d*', src_ref),
            '\'{}\',{},{}'.format(self.gcs.spheroidName,
                                  spheroid['Semimajor_Axis'],
                                  spheroid['Inverse_Flattening'])))
        revision.append(('GCS_North_American_1983', self.gcs.name))
        revision.append(('D_North_American_1983', self.gcs.datumName))
        for i in revision:
            src_ref = src_ref.replace(*i)
        new_ref = arcpy.SpatialReference()
        new_ref.loadFromString(src_ref)
        return new_ref


    def __modify_secant_parameters(self, spatial_ref, extent):

        """
        Modifies standard parallels and central meridan for spatial
        references.
        """
        src_ref = spatial_ref.exportToString()
        revision = list()
        secant = self.__get_secant_parameters(extent)
        for i, j in secant.items():
            revision.append(
                (self.__find('\'{}\',[-+]?\d*[.]\d*'.format(i), src_ref),
                 '\'{}\',{}'.format(i, j)))
        for i in revision:
            src_ref = src_ref.replace(*i)
        new_ref = arcpy.SpatialReference()
        new_ref.loadFromString(src_ref)
        return new_ref


    def __find(self, pattern, target):

        """
        Returns the first result of a REGEX pattern against a target string.
        """
        return re.compile(pattern).findall(target)[0]
    

    def __get_main_extent(self):

        """
        Returns the extent of a feature class, projected to the specified GCS
        for use in projection parameters.
        """
        polygon = arcpy.Describe(self.features).extent.polygon
        project_polygon = polygon.projectAs(self.gcs, self.main_transformation)
        return project_polygon.extent


    def __get_geometry_extent(self, geometry):

        """
        Returns the extent of a geometry object, projected to a specified GCS
        for use in projection parameters.
        """
        polygon = geometry.extent.polygon
        project_polygon = polygon.projectAs(self.gcs, self.main_transformation)
        return project_polygon.extent


    def __get_main_transformation(self):

        """
        Returns the primary transformation covering the extent of the main
        feature class.
        """
        feature_info = arcpy.Describe(self.features)
        parameters = {'from_sr': feature_info.spatialReference,
                      'to_sr': self.gcs,
                      'extent': feature_info.extent}
        return self.__get_transformation(**parameters)


    def __get_transformation(self, from_sr, to_sr, extent=None):

        """
        Returns the primary transformation covering a specified extent.
        """
        parameters = {'from_sr': from_sr, 'to_sr': to_sr, 'extent': extent}
        transformations = arcpy.ListTransformations(**parameters)
        if transformations:
            return transformations[0]
        else:
            return None


    def __get_secant_parameters(self, extent):

        """
        Returns standard parallels and central meridian based on specified
        extent (requires extent with GCS rather than PCS).
        """
        x_min, x_max = extent.XMin, extent.XMax
        y_min, y_max = extent.YMin, extent.YMax
        parallel_coeff = (y_max - y_min) / 6.0
        return {'Central_Meridian': ((x_max - x_min) / 2.0) + x_min,
                'Standard_Parallel_1': y_min + parallel_coeff,
                'Standard_Parallel_2': y_max - parallel_coeff}


    def __get_spheroid_parameters(self):

        """
        Returns the spheroid parameters for a projection based on the main GCS
        of the object.
        """
        semiMajAxis = self.gcs.semiMajorAxis
        semiMinAxis = self.gcs.semiMinorAxis
        inverse_flattening = semiMajAxis / (semiMajAxis - semiMinAxis)
        return {'Spheroid': self.gcs.spheroidName,
                'Semimajor_Axis': semiMajAxis,
                'Inverse_Flattening': inverse_flattening}
    

    def __msg(self, text):

        """
        Provides dual messaging in both Python console and geoprocessing
        messages.
        """
        print(text)
        arcpy.AddMessage(text)
        return


    def __calc_sq_m_to_acres(self, value):

        """
        Returns acreage calculation based on square meters to the highest
        accuracy and precision available.
        """
        return value / (4046.0 + (13525426.0 / 15499969.0))


    def __create_acreage_field(self):

        """
        Creates a field within the feature class to contain acreage
        information.  If the field already exists, it will be deleted and
        recreated.
        """
        field_list = [i.name for i in arcpy.ListFields(
            self.features, 'GIS_Acres')]
        if field_list:
            arcpy.DeleteField_management(self.features, 'GIS_Acres')
        parameters = {'in_table': self.features, 'field_name': 'GIS_Acres',
                      'field_alias': 'GIS Acres', 'field_type': 'DOUBLE'}
        arcpy.AddField_management(**parameters)
        return


    def __calculate_acreage(self):

        """
        Calculates the acreage of features in a feature class. Using an
        UpdateCursor, retrieves feature geometry and instantiates projected
        Geometry objects to derive area, then calculates acreage. Updates the
        feature row with the calculated acreage value.
        """
        self.__msg('Calculating acreage...')
        self.__create_acreage_field()
        fields = ['SHAPE@', 'GIS_Acres']
        with arcpy.da.UpdateCursor(self.features, fields) as cursor:
            if not self.local:
                for row in cursor:
                    geometry = row[0].projectAs(self.main_spatial_ref,
                                                self.main_transformation)
                    row[1] = self.__calc_sq_m_to_acres(geometry.area)
                    self.__msg('- {} acres'.format(row[1]))
                    cursor.updateRow(row)
            else:
                for row in cursor:
                    src_geom = row[0]
                    src_extent = self.__get_geometry_extent(src_geom)
                    src_ref = self.__modify_secant_parameters(
                        self.main_spatial_ref, src_extent)
                    geometry = src_geom.projectAs(
                        src_ref, self.main_transformation)
                    row[1] = self.__calc_sq_m_to_acres(geometry.area)
                    self.__msg('{} acres'.format(row[1]))
                    cursor.updateRow(row)
        del cursor # Necessary even with exit conditions in with statement
        return


if __name__ == '__main__':

    try:
        received_params = list(arcpy.GetParameter(i) for i in
                               range(arcpy.GetArgumentCount()))
        input_params = {'feature_class': str(received_params[0]),
                        'local': received_params[1]}
        calc_object = AcreageCalculator(**input_params)
    except:
        print('Unable to receive geoprocessing parameters from ArcGIS')  

