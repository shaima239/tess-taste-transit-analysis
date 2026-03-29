import numpy as np
import matplotlib.pyplot as plt
import pickle
from astropy.io import fits
from astropy import coordinates as coord, units as u
from astropy.time import Time

class AperturePhotometry:
    data_dir = './group15_WASP-12_20191229/'

    def __init__(self):
        """
        Initializing the class by loading constants,
        calibration frames, and data configurations.
        """        

        self.data_path = self.data_dir

        # Calibration and instrument parameters. 
        self.readout_noise = 7.1  # [e-] photoelectrons
        self.gain = 1.91  # [e-/ADU]
        self.bias_std = 1.3  # [e-] photoelectrons

        # Loading median bias and its errors.
        self.median_bias = pickle.load(open('median_bias.p', 'rb'))
        self.median_bias_errors = pickle.load(open('median_bias_error.p', 'rb'))

        # Loading median flat and its error.
        self.median_normalized_flat = pickle.load(open('median_normalized_flat.p', 'rb'))
        self.median_normalized_flat_errors = pickle.load(open('median_normalized_flat_errors.p', 'rb'))

        # Loading cience file list
        self.science_path = self.data_path + 'science/'
        self.science_list = np.genfromtxt(self.science_path + 'science.list', dtype=str)
        self.science_size = len(self.science_list)

        # Generating a meshgrid for pixel coordinates.
        ylen, xlen = np.shape(self.median_bias)
        X_axis = np.arange(0, xlen, 1)
        Y_axis = np.arange(0, ylen, 1)
        self.X, self.Y = np.meshgrid(X_axis, Y_axis)
        self.X_axis = X_axis
        self.Y_axis = Y_axis

        # Defining our target and the observatory's location.
        self.target = coord.SkyCoord("06:30:32.79", "+29:40:20", unit=(u.hourangle, u.deg), frame='icrs')
        self.observatory_location = ('45.872d', '11.569d')  # Asiago observatory coordinates

    def provide_aperture_parameters(self, sky_inner_radius, sky_outer_radius, aperture_radius, x_initial, y_initial):
        """
        Providing parameters for aperture photometry and the annulus around our 
        star and the reference stars.
        
        Arguments x_initial and y_initial: Initial coordinates of the star [pixels].
        """
        
        self.sky_inner_radius = sky_inner_radius
        self.sky_outer_radius = sky_outer_radius
        self.aperture_radius = aperture_radius
        self.x_initial = x_initial
        self.y_initial = y_initial

    def correct_science_frame(self, science_frame):
        """
        Correcting the science frame for bias and flat-field.

        Returns:
            Bias and flat-corrected science image and its errors
        """        
        science_debiased = science_frame - self.median_bias
        science_corrected = science_debiased / self.median_normalized_flat

        # Coputing errors for the corrected frame 
        science_debiased_errors = np.sqrt(self.readout_noise**2 + self.bias_std**2 + science_debiased)

        # Generating a mask to ensure the values in the pixle are not zero.
        valid = (science_debiased != 0) & (self.median_normalized_flat != 0)
        
        science_corrected_errors = np.zeros_like(science_corrected)
        science_corrected_errors[valid] = science_corrected[valid] * np.sqrt(
            (science_debiased_errors[valid] / science_debiased[valid])**2 +
            (self.median_normalized_flat_errors[valid] / self.median_normalized_flat[valid])**2
        )
        science_corrected_errors[~valid] = 0.0

        return science_corrected, science_corrected_errors

    def compute_centroid(self, science_frame, x_target_initial, y_target_initial, maximum_number_of_iterations=20):
        """
        Refining the centroid's position using an iterative method.

        Returns:
            x_refined and y_refined: Refined coordinates of the centroid.
        """
        
        for _ in range(maximum_number_of_iterations):
            target_distance = np.sqrt((self.X - x_target_initial)**2 + (self.Y - y_target_initial)**2)
            annulus_sel = (target_distance < self.sky_inner_radius)

            total_flux = np.sum(science_frame[annulus_sel])
            if total_flux == 0:
                break

            weighted_X = np.sum(science_frame[annulus_sel] * self.X[annulus_sel])
            weighted_Y = np.sum(science_frame[annulus_sel] * self.Y[annulus_sel])

            x_refined = weighted_X / total_flux
            y_refined = weighted_Y / total_flux

            if np.abs(x_refined - x_target_initial) < 0.1 and np.abs(y_refined - y_target_initial) < 0.1:
                break

            x_target_initial, y_target_initial = x_refined, y_refined

        return x_refined, y_refined

    def compute_sky_background(self, science_frame, science_frame_errors, x_pos, y_pos):
        """
        Computing the sky background and its error.

        Returns:
            Median sky background flux and its errors.
        """
        target_distance = np.sqrt((self.X - x_pos)**2 + (self.Y - y_pos)**2)
        annulus_selection = (target_distance > self.sky_inner_radius) & (target_distance <= self.sky_outer_radius)

        sky_flux_median = np.median(science_frame[annulus_selection])
        N = np.sum(annulus_selection)
        sky_flux_error = np.sqrt(np.sum(science_frame_errors[annulus_selection] ** 2)) / N

        return sky_flux_median, sky_flux_error

    def determine_FWHM_axis(self, reference_axis, normalized_cumulative_distribution):
        """
        Determining the Full Width at Half Maximum (FWHM) along a specific axis.
    
        Args:
            reference_axis: Pixel positions along the axis (e.g., X or Y).
            normalized_cumulative_distribution: Normalized cumulative sum of pixel intensities along the axis.
    
        Returns:
            FWHM: Full Width at Half Maximum (FWHM) in pixels.
        """
        # Finding the closest points to NCD=0.15865 (-1σ) and NCD=0.84135 (+1σ)
        NCD_index_left = np.argmin(np.abs(normalized_cumulative_distribution - 0.15865))
        NCD_index_right = np.argmin(np.abs(normalized_cumulative_distribution - 0.84135))
    
        # Fiting polynomial to refine positions
        p_left = np.polynomial.Polynomial.fit(
            normalized_cumulative_distribution[NCD_index_left - 1: NCD_index_left + 2],
            reference_axis[NCD_index_left - 1: NCD_index_left + 2],
            deg=2
        )
        pixel_left = p_left(0.15865)
    
        p_right = np.polynomial.Polynomial.fit(
            normalized_cumulative_distribution[NCD_index_right - 1: NCD_index_right + 2],
            reference_axis[NCD_index_right - 1: NCD_index_right + 2],
            deg=2
        )
        pixel_right = p_right(0.84135)
    
        # Converting to FWHM
        FWHM_factor = 2 * np.sqrt(2 * np.log(2))  # = 2.35482
        FWHM = (pixel_right - pixel_left) / 2. * FWHM_factor
    
        return FWHM


    def compute_fwhm(self, science_frame, x_centroid, y_centroid, fit_radius=10):
        """
        Computing the FWHM of the star in the science frame.
        
        Args:
            science_frame: 2D array of the science image data.
            x_centroid: x-coordinate of the centroid.
            y_centroid: y-coordinate of the centroid.
            fit_radius: Radius around the centroid for fitting.
        
            Returns:
                FWHM_x: FWHM in the x-direction.
                FWHM_y: FWHM in the y-direction.
        """
        # Extracting slices through the centroid.
        x_slice = science_frame[int(y_centroid), int(x_centroid - fit_radius):int(x_centroid + fit_radius)]
        y_slice = science_frame[int(y_centroid - fit_radius):int(y_centroid + fit_radius), int(x_centroid)]
        
        # Computing cumulative sums.
        cumulative_sum_x = np.cumsum(x_slice) / np.sum(x_slice)
        cumulative_sum_y = np.cumsum(y_slice) / np.sum(y_slice)
        
        # Creating reference axes.
        x_range = np.arange(len(x_slice))
        y_range = np.arange(len(y_slice))
        
        # Calculating FWHM using the normalized cumulative distribution
        FWHM_x = self.determine_FWHM_axis(x_range, cumulative_sum_x)
        FWHM_y = self.determine_FWHM_axis(y_range, cumulative_sum_y)
        
        return FWHM_x, FWHM_y



    def aperture_photometry(self):
        """
        Performing aperture photometry for all science images in the dataset.
        """
        self.airmass = np.empty(self.science_size)
        self.exptime = np.empty(self.science_size)
        self.julian_date = np.empty(self.science_size)

        self.aperture = np.empty(self.science_size)
        self.aperture_errors = np.empty(self.science_size)
        self.sky_background = np.empty(self.science_size)
        self.sky_background_errors = np.empty(self.science_size)

        self.x_position = np.empty(self.science_size)
        self.y_position = np.empty(self.science_size)
        self.x_fwhm = np.empty(self.science_size)
        self.y_fwhm = np.empty(self.science_size)

        for i, science_name in enumerate(self.science_list):
            try:
                science_fits = fits.open(self.science_path + science_name)
                hdr = science_fits[0].header
                science_data = science_fits[0].data * self.gain
                science_fits.close()

                self.airmass[i] = hdr['AIRMASS']
                self.exptime[i] = hdr['EXPTIME']
                self.julian_date[i] = hdr['JD']

                # Correcting science frame
                science_corrected, science_corrected_errors = self.correct_science_frame(science_data)

                # Computing centroid
                x_refined, y_refined = self.compute_centroid(science_corrected, self.x_initial, self.y_initial)

                # Computing sky background
                sky_flux_median, sky_flux_error = self.compute_sky_background(science_corrected, science_corrected_errors, x_refined, y_refined)
                self.sky_background[i] = sky_flux_median
                self.sky_background_errors[i] = sky_flux_error

                # Subtracting sky and computing aperture photometry
                science_sky_corrected = science_corrected - sky_flux_median
                target_distance = np.sqrt((self.X - x_refined)**2 + (self.Y - y_refined)**2)
                aperture_selection = (target_distance < self.aperture_radius)

                num_pixels_in_aperture = np.sum(aperture_selection)

                self.aperture[i] = np.sum(science_sky_corrected[aperture_selection])
                photon_noise = np.sum(science_sky_corrected[aperture_selection])
                readout_noise = self.readout_noise * np.sqrt(num_pixels_in_aperture)
                self.aperture_errors[i] = np.sqrt(photon_noise + sky_flux_error**2 * num_pixels_in_aperture + readout_noise**2)

                self.x_position[i] = x_refined
                self.y_position[i] = y_refined

                # Compute FWHM
                fwhm_x, fwhm_y = self.compute_fwhm(science_sky_corrected, x_refined, y_refined, self.sky_inner_radius)
                self.x_fwhm[i] = fwhm_x
                self.y_fwhm[i] = fwhm_y

            except Exception as e:
                print(f"Error processing {science_name}: {e}")
                self.aperture[i] = np.nan
                self.aperture_errors[i] = np.nan

        # Converting JD to BJD_TDB
        tm = Time(self.julian_date, format='jd', scale='utc', location=('45.872d', '11.569d')) #self.julian_date[i_science]
        llt_bary = tm.light_travel_time(self.target)
        self.bjd_tdb = (tm.tdb + llt_bary).to_value('jd')

