import { useI18nStore, translations } from '../lib/i18n';
import { LanguageSwitcher } from './LanguageSwitcher';

export function Header() {
  const { language } = useI18nStore();
  const t = translations[language];

  return (
    <header className="sticky top-0 z-50 w-full bg-gradient-to-r from-blue-900 to-indigo-800 shadow-md">
      <div className="container mx-auto px-4 py-3 flex items-center justify-between">
        <div className="flex items-center space-x-2">
          <svg 
            className="w-8 h-8 text-blue-200" 
            viewBox="0 0 24 24" 
            fill="none" 
            stroke="currentColor" 
            strokeWidth="2" 
            strokeLinecap="round" 
            strokeLinejoin="round"
          >
            <path d="M12 2L2 7l10 5 10-5-10-5z" />
            <path d="M2 17l10 5 10-5" />
            <path d="M2 12l10 5 10-5" />
          </svg>
          <span className="text-xl font-bold text-white">
            {language === 'zh' ? '木薯泛基因组' : 'Cassava Pangenome'}
          </span>
        </div>
        
        <nav className="hidden md:flex items-center space-x-6">
          <a href="#" className="text-white hover:text-blue-200 transition-colors">
            {t.nav.home}
          </a>
          <a href="#" className="text-white hover:text-blue-200 transition-colors">
            {t.nav.github}
          </a>
          <a href="#" className="text-white hover:text-blue-200 transition-colors">
            {t.nav.documentation}
          </a>
        </nav>
        
        <div className="flex items-center space-x-4">
          <LanguageSwitcher />
          <button className="md:hidden text-white">
            <svg 
              className="w-6 h-6" 
              fill="none" 
              stroke="currentColor" 
              viewBox="0 0 24 24"
            >
              <path 
                strokeLinecap="round" 
                strokeLinejoin="round" 
                strokeWidth="2" 
                d="M4 6h16M4 12h16M4 18h16"
              />
            </svg>
          </button>
        </div>
      </div>
    </header>
  );
}
