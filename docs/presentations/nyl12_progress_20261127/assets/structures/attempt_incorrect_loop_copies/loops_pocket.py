"""Render the actual A/E pocket with distinct candidate loop regions."""

from pathlib import Path
import pymol2


def main() -> None:
    """Render an isolated mature protomer without changing its coordinates."""
    here = Path(__file__).resolve().parent
    with pymol2.PyMOL() as instance:
        cmd = instance.cmd
        cmd.load(str(here / "assembly_loops.pse"))
        cmd.set("max_threads", 2)
        cmd.hide("everything", "mol and not chain A+E+I")
        cmd.hide(
            "cartoon",
            "mol and ((chain A and resi 95-103) or (chain E and resi 289-300))",
        )
        for obj, selection, color in [
            ("loop_alpha", "chain A and resi 95-103", "loop_alpha"),
            ("loop_beta", "chain E and resi 289-300", "loop_beta"),
        ]:
            cmd.create(obj, "mol and (" + selection + ")")
            cmd.hide("everything", obj)
            cmd.show("cartoon", obj)
            cmd.cartoon("tube", obj)
            cmd.set("cartoon_tube_radius", 0.5, obj)
            cmd.color(color, obj)
        cmd.zoom("mol and chain A+E+I", buffer=3)
        cmd.png(str(here / "loops_pocket.png"), 1800, 1600, ray=1, quiet=1)
        cmd.save(str(here / "loops_pocket.pse"))


if __name__ == "__main__":
    main()
