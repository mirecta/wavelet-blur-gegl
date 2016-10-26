/* This file is an image processing operation for GEGL
 * GEGL is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * GEGL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GEGL; if not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2016 Miroslav Talasek <miroslav.talasek@seznam.cz>
 *
 * Wavelet blur used in wavelet decompose filter
 *  theory is from original wavelet plugin
 *
 */

#include "config.h"
#include <glib/gi18n-lib.h>

#ifdef GEGL_PROPERTIES


property_double (radius, _("Radius"), 1.0)
    description (_("Radius of wavelet blur"))
    value_range (0.0, 1500.0)
    ui_range    (0.0, 256.0)
    ui_gamma    (3.0)
    ui_meta     ("unit", "pixel-distance")
    ui_meta     ("radius", "blur")


#else

#define GEGL_OP_AREA_FILTER
#define GEGL_OP_NAME     wavelet_blur
#define GEGL_OP_C_SOURCE wavelet-blur.c

#include "gegl-op.h"
#include <math.h>
#include <stdio.h>


static gint
wav_gen_convolve_matrix (gdouble   sigma,
                         gdouble **cmatrix_p);



static gint
wav_calc_convolve_matrix_length (gdouble rad)
{
  return ceil (rad) * 2 + 1;
}

static gint
wav_gen_convolve_matrix (gdouble   sigma,
                         gdouble **cmatrix_p)
{
  gint     matrix_length;
  gdouble *cmatrix;

  matrix_length = wav_calc_convolve_matrix_length (sigma);
  cmatrix = g_new (gdouble, matrix_length);
  if (!cmatrix)
    return 0;

  if (matrix_length == 1)
    {
      cmatrix[0] = 1;
    }
  else
    {
      gint i,x;
      gdouble sum = 0;

      for (i=0; i<matrix_length; i++)
        {
          if (i == 0 || i == matrix_length - 1)
            {
              cmatrix[i] = 0.25;
            }
          else if (i == matrix_length /2)
            {
              cmatrix[i] = 0.5;
            }
          else
            {
              cmatrix[i] = 0;
            }
        }
    }

  *cmatrix_p = cmatrix;
  return matrix_length;
}

static inline void
wav_get_mean_pixel_1D (gfloat  *src,
                       gfloat  *dst,
                       gint     components,
                       gdouble *cmatrix,
                       gint     matrix_length)
{
  gint    i, c;
  gint    offset;
  gdouble acc[components];

  for (c = 0; c < components; ++c)
    acc[c] = 0;

  offset = 0;

  for (i = 0; i < matrix_length; i++)
    {
      for (c = 0; c < components; ++c)
        acc[c] += src[offset++] * cmatrix[i];
    }

  for (c = 0; c < components; ++c)
    dst[c] = acc[c];
}

static void 
wav_fix_bordes (GeglBuffer *buf,
                gint buf_length,
                gint radius)
{
  gint     i, c;
  
  //left side
  for (i = 0; i < radius; ++i)
    {
     for (c = 0; c < components; ++c)
       {
        buf[i*components + c] = buf[2*radius*components - i*components + c];
       }
    }
  
}
                

static void
wav_hor_blur (GeglBuffer          *src,
              GeglBuffer          *dst,
              const GeglRectangle *dst_rect,
              gdouble             *cmatrix,
              gint                 matrix_length)
{
  gint        u, v;
  const gint  radius = matrix_length / 2;
  const Babl *format = babl_format ("RaGaBaA float");

  GeglRectangle write_rect = {dst_rect->x, dst_rect->y, dst_rect->width, 1};
  gfloat *dst_buf     = gegl_malloc (write_rect.width * sizeof(gfloat) * 4);

  GeglRectangle read_rect = {dst_rect->x - radius, dst_rect->y, dst_rect->width + radius, 1};
  gfloat *src_buf    = gegl_malloc (read_rect.width * sizeof(gfloat) * 4);

  for (v = 0; v < dst_rect->height; v++)
    {
      gint offset     = 0;
      read_rect.y     = dst_rect->y + v;
      write_rect.y    = dst_rect->y + v;
      gegl_buffer_get (src, &read_rect, 1.0, format, src_buf, GEGL_AUTO_ROWSTRIDE, GEGL_ABYSS_NONE);
      //fix edges
      for (u = 0; u < dst_rect->width; u++)
        {
          wav_get_mean_pixel_1D (src_buf + offset,
                                 dst_buf + offset,
                                 4,
                                 cmatrix,
                                 matrix_length);
          offset += 4;
        }

      gegl_buffer_set (dst, &write_rect, 0, format, dst_buf, GEGL_AUTO_ROWSTRIDE);
    }

  gegl_free (src_buf);
  gegl_free (dst_buf);
}

static void
wav_ver_blur (GeglBuffer          *src,
              GeglBuffer          *dst,
              const GeglRectangle *dst_rect,
              gdouble             *cmatrix,
              gint                 matrix_length)
{
  gint        u,v;
  const gint  radius = matrix_length / 2;
  const Babl *format = babl_format ("RaGaBaA float");

  GeglRectangle write_rect = {dst_rect->x, dst_rect->y, 1, dst_rect->height};
  gfloat *dst_buf    = gegl_malloc (write_rect.height * sizeof(gfloat) * 4);

  GeglRectangle read_rect  = {dst_rect->x, dst_rect->y - radius, 1, dst_rect->height + radius};
  gfloat *src_buf    = gegl_malloc (read_rect.height * sizeof(gfloat) * 4);

  for (u = 0; u < dst_rect->width; u++)
    {
      gint offset     = 0;
      read_rect.x     = dst_rect->x + u;
      write_rect.x    = dst_rect->x + u;
      gegl_buffer_get (src, &read_rect, 1.0, format, src_buf, GEGL_AUTO_ROWSTRIDE, GEGL_ABYSS_NONE);
      //fix edges
      for (v = 0; v < dst_rect->height; v++)
        {
          wav_get_mean_pixel_1D (src_buf + offset,
                                 dst_buf + offset,
                                 4,
                                 cmatrix,
                                 matrix_length);
          offset += 4;
        }

      gegl_buffer_set (dst, &write_rect, 0, format, dst_buf, GEGL_AUTO_ROWSTRIDE);
    }

  gegl_free (src_buf);
  gegl_free (dst_buf);
}

static void
prepare (GeglOperation *operation)
{
  GeglOperationAreaFilter *area = GEGL_OPERATION_AREA_FILTER (operation);
  GeglProperties          *o    = GEGL_PROPERTIES (operation);

  gfloat fir_radius_x = fir_calc_convolve_matrix_length (o->std_dev_x) / 2;
  gfloat fir_radius_y = fir_calc_convolve_matrix_length (o->std_dev_y) / 2;
  gfloat iir_radius_x = o->std_dev_x * RADIUS_SCALE;
  gfloat iir_radius_y = o->std_dev_y * RADIUS_SCALE;

  /* XXX: these should be calculated exactly considering o->filter, but we just
   * make sure there is enough space */
  area->left = area->right = ceil (MAX (fir_radius_x, iir_radius_x));
  area->top = area->bottom = ceil (MAX (fir_radius_y, iir_radius_y));

  gegl_operation_set_format (operation, "input",
                             babl_format ("RaGaBaA float"));
  gegl_operation_set_format (operation, "output",
                             babl_format ("RaGaBaA float"));
}

static gboolean
process (GeglOperation       *operation,
         GeglBuffer          *input,
         GeglBuffer          *output,
         const GeglRectangle *result,
         gint                 level)
{
  GeglRectangle rect;
  GeglBuffer *temp;
  GeglOperationAreaFilter *op_area = GEGL_OPERATION_AREA_FILTER (operation);
  GeglProperties          *o       = GEGL_PROPERTIES (operation);

  GeglRectangle temp_extend;
  gdouble       B, b[4];
  gdouble      *cmatrix;
  gint          cmatrix_len;
  gboolean      horizontal_irr;
  gboolean      vertical_irr;

  rect.x      = result->x - op_area->left;
  rect.width  = result->width + op_area->left + op_area->right;
  rect.y      = result->y - op_area->top;
  rect.height = result->height + op_area->top + op_area->bottom;

  if (o->filter == GEGL_GAUSSIAN_BLUR_FILTER_IIR)
    {
      horizontal_irr = TRUE;
      vertical_irr   = TRUE;
    }
  else if (o->filter == GEGL_GAUSSIAN_BLUR_FILTER_FIR)
    {
      horizontal_irr = FALSE;
      vertical_irr   = FALSE;
    }
  else /* GEGL_GAUSSIAN_BLUR_FILTER_AUTO */
    {
      horizontal_irr = o->std_dev_x > 1.0;
      vertical_irr   = o->std_dev_y > 1.0;
    }

  if (gegl_operation_use_opencl (operation) && !(horizontal_irr | vertical_irr))
    if (cl_process(operation, input, output, result))
      return TRUE;

  gegl_rectangle_intersect (&temp_extend, &rect, gegl_buffer_get_extent (input));
  temp_extend.x      = result->x;
  temp_extend.width  = result->width;
  temp = gegl_buffer_new (&temp_extend, babl_format ("RaGaBaA float"));

  if (horizontal_irr)
    {
      iir_young_find_constants (o->std_dev_x, &B, b);
      iir_young_hor_blur (input, &rect, temp, &temp_extend, B, b);
    }
  else
    {
      cmatrix_len = fir_gen_convolve_matrix (o->std_dev_x, &cmatrix);
      fir_hor_blur (input, temp, &temp_extend, cmatrix, cmatrix_len);
      g_free (cmatrix);
    }

  if (vertical_irr)
    {
      iir_young_find_constants (o->std_dev_y, &B, b);
      iir_young_ver_blur (temp, &rect, output, result, B, b);
    }
  else
    {
      cmatrix_len = fir_gen_convolve_matrix (o->std_dev_y, &cmatrix);
      fir_ver_blur (temp, output, result, cmatrix, cmatrix_len);
      g_free (cmatrix);
    }

  g_object_unref (temp);
  return  TRUE;
}


static void
gegl_op_class_init (GeglOpClass *klass)
{
  GeglOperationClass       *operation_class;
  GeglOperationFilterClass *filter_class;

  operation_class = GEGL_OPERATION_CLASS (klass);
  filter_class    = GEGL_OPERATION_FILTER_CLASS (klass);

  operation_class->prepare        = prepare;
  operation_class->opencl_support = TRUE;

  filter_class->process           = process;

  gegl_operation_class_set_keys (operation_class,
    "name",        "gegl:gaussian-blur-old",
    "title",       _("Gaussian Blur"),
    "categories",  "blur",
    "description", _("Each result pixel is the average of the neighbouring pixels weighted by a normal distribution with specified standard deviation."),
    NULL);
}

#endif
