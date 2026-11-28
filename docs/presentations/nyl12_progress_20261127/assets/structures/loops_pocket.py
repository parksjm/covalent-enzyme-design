"""Render actual neighboring subunit loops around source-bound L4."""

from pathlib import Path
import pymol2


def main() -> None:
    """Render full assembly and a local view without changing coordinates."""
    here = Path(__file__).resolve().parent
    with pymol2.PyMOL() as instance:
        cmd = instance.cmd
        cmd.load(str(here / "assembly_loops.pse"))
        cmd.set("max_threads", 2)
        for chains, color in [
            ("A+E", "green"),
            ("B+F", "cyan"),
            ("C+G", "magenta"),
            ("D+H", "yellow"),
        ]:
            cmd.color(color, "mol and chain " + chains)
        cmd.color("substrate_c", "mol and chain I and elem C")
        cmd.color("red", "mol and (chain I or (chain E and resi 233)) and elem O")
        cmd.color("blue", "mol and (chain I or (chain E and resi 233)) and elem N")
        for obj, selection, color in [
            ("loop_alpha", "chain D and resi 95-103", "loop_alpha"),
            ("loop_beta", "chain F and resi 289-300", "loop_beta"),
        ]:
            cmd.hide("cartoon", "mol and (" + selection + ")")
            cmd.create(obj, "mol and (" + selection + ")")
            cmd.hide("everything", obj)
            cmd.show("cartoon", obj)
            cmd.cartoon("tube", obj)
            cmd.set("cartoon_tube_radius", 0.5, obj)
            cmd.color(color, obj)
        cmd.zoom("mol", buffer=3)
        cmd.png(str(here / "assembly_loops.png"), 1800, 1600, ray=1, quiet=1)
        cmd.save(str(here / "assembly_loops.pse"))
        cmd.set("cartoon_transparency", 0.55, "mol")
        cmd.zoom("loop_alpha or loop_beta or (mol and chain I)", buffer=4)
        cmd.png(str(here / "loops_pocket.png"), 1800, 1600, ray=1, quiet=1)
        cmd.save(str(here / "loops_pocket.pse"))


if __name__ == "__main__":
    main()
