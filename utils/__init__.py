import bpy


def update_screen():
    try:
        for area in bpy.context.screen.areas:
            area.tag_redraw()
    except BaseException:
        ...
