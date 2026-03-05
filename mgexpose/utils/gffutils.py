""" Module docstring """
def get_source_column(source_db=None,):
    """ docstring """
    return ("proMGE", f"proMGE_{source_db}")[bool(source_db)]


def get_attrib_str(attribs):
    """ docstring """
    return ";".join(f"{item[0]}={item[1]}" for item in attribs.items() if item[1])


def parse_gff_attribs(attrib_str):
    """ docstring """
    try:
        attribs = dict(item.split("=") for item in attrib_str.split(";"))
    except Exception as exc:
        raise ValueError(f"problem parsing gff attribute string: {attrib_str}") from exc
    return attribs
