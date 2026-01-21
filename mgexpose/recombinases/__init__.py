""" Module for recombinase functions. """

from collections import Counter

from .recombinases import MGE_ALIASES, MgeRule


def parse_recombinase_string(recombinase_str):
	try:
		recombinases = Counter(
			dict(
				(key, int(value))
				for key, value in (
					item.split(":") for item in recombinase_str.split(",")
				)
			)
		)
	except Exception as exc:
		raise ValueError(f"recombinase string weird? {recombinase_str}") from exc
	
	return recombinases
