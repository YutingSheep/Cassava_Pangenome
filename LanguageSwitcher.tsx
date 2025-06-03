import { useI18nStore, type Language } from '../lib/i18n';

export function LanguageSwitcher() {
  const { language, setLanguage } = useI18nStore();

  const toggleLanguage = () => {
    setLanguage(language === 'zh' ? 'en' : 'zh');
  };

  return (
    <button
      onClick={toggleLanguage}
      className="flex items-center justify-center px-3 py-1 text-sm font-medium rounded-md bg-white/10 hover:bg-white/20 text-white transition-colors"
      aria-label={language === 'zh' ? "Switch to English" : "切换到中文"}
    >
      {language === 'zh' ? 'EN' : '中文'}
    </button>
  );
}
