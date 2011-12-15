! -*- f77 -*-
      subroutine get_elem_config(out_nsm, out_ncm)

      implicit none

      integer out_nsm, out_ncm

      include "elem_config.inc"

      out_nsm = NELMAX+2
      out_ncm = NIONMAX

      return

      end
