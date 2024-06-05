from pathlib import Path
import importlib


def compile_translation(translations: tuple[tuple[str, str]]):
    t = {}
    for item in translations:
        if len(item) < 2:
            continue
        context = None if len(item) == 2 else item[2]
        source, translation = item[:2]
        t[(context, source)] = translation
    return t


def load_translations():
    translations_dir = Path(__file__).parent / "translations"
    translations_dict = {}

    for translation_file in translations_dir.glob("*.py"):
        if translation_file.stem != "__init__":
            language_code = translation_file.stem
            locale = language_code.replace('-', '_')
            translation_module = importlib.import_module(f".translations.{language_code}", package=__package__)
            if not hasattr(translation_module, "translations"):
                continue
            translations = getattr(translation_module, "translations")
            translations_dict[locale] = compile_translation(translations)

    return translations_dict
