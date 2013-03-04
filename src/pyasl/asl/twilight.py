
def twilightName(sunAlt):
  """
    Determine the type of twilight for a given altitude of the Sun.
    
    Definitions of twilight:
      - Solar altitude > 0    : day
      - Solar altitude <= 0   : civil twilight
      - Solar altitude <= -6  : nautical twilight
      - Solar altitude <= -12 : astronomical twilight
      - Solar altitude <= -18 : night

    Parameters
    ----------
    sunAlt : float
        Altitude of the Sun in degrees.
         
    Returns
    -------
    Twilight flag : string, {"day", "civil twilight", "nautical twilight", "astron. twilight", "night"}
        Specifies the type of twilight corresponding to the
        given solar altitude.
         
        The twilight flag is either "day", "civil twilight", "nautical twilight",
        "astron. twilight", or "night".
  """  
  
  if sunAlt <= -18.0: return "night"
  if sunAlt <= -12.0: return "astron. twilight"
  if sunAlt <=  -6.0: return "nautical twilight"
  if sunAlt <=   0.0: return "civil twilight"
  return "day"