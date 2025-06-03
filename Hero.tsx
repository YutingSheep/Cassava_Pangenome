import { useI18nStore, translations } from '../lib/i18n';

export function Hero() {
  const { language } = useI18nStore();
  const t = translations[language];

  return (
    <section className="relative bg-gradient-to-b from-blue-900 to-indigo-950 text-white py-16 md:py-24 overflow-hidden">
      {/* 装饰性DNA双螺旋背景 */}
      <div className="absolute inset-0 opacity-10">
        <svg width="100%" height="100%" viewBox="0 0 100 100" preserveAspectRatio="none">
          <defs>
            <pattern id="dna-pattern" patternUnits="userSpaceOnUse" width="100" height="100" patternTransform="rotate(45)">
              <path d="M0,0 L100,0" stroke="white" strokeWidth="1" fill="none" strokeDasharray="5,10" />
              <path d="M0,20 L100,20" stroke="white" strokeWidth="1" fill="none" strokeDasharray="10,5" />
              <path d="M0,40 L100,40" stroke="white" strokeWidth="1" fill="none" strokeDasharray="5,10" />
              <path d="M0,60 L100,60" stroke="white" strokeWidth="1" fill="none" strokeDasharray="10,5" />
              <path d="M0,80 L100,80" stroke="white" strokeWidth="1" fill="none" strokeDasharray="5,10" />
            </pattern>
          </defs>
          <rect width="100%" height="100%" fill="url(#dna-pattern)" />
        </svg>
      </div>
      
      <div className="container mx-auto px-4 relative z-10">
        <div className="flex flex-col md:flex-row items-center">
          <div className="md:w-1/2 mb-8 md:mb-0">
            <h1 className="text-4xl md:text-5xl font-bold mb-4 bg-clip-text text-transparent bg-gradient-to-r from-blue-200 to-indigo-100">
              {t.hero.title}
            </h1>
            <p className="text-xl text-blue-200 mb-6">{t.hero.subtitle}</p>
            <p className="text-gray-300 mb-8 max-w-lg">{t.hero.description}</p>
            <div className="flex flex-wrap gap-4">
              <a 
                href="#modules" 
                className="px-6 py-3 bg-blue-600 hover:bg-blue-700 text-white font-medium rounded-lg transition-colors"
              >
                {t.hero.getStarted}
              </a>
              <a 
                href="https://github.com" 
                target="_blank" 
                rel="noopener noreferrer"
                className="px-6 py-3 bg-gray-700 hover:bg-gray-800 text-white font-medium rounded-lg transition-colors flex items-center gap-2"
              >
                <svg className="w-5 h-5" viewBox="0 0 24 24" fill="currentColor">
                  <path d="M12 0c-6.626 0-12 5.373-12 12 0 5.302 3.438 9.8 8.207 11.387.599.111.793-.261.793-.577v-2.234c-3.338.726-4.033-1.416-4.033-1.416-.546-1.387-1.333-1.756-1.333-1.756-1.089-.745.083-.729.083-.729 1.205.084 1.839 1.237 1.839 1.237 1.07 1.834 2.807 1.304 3.492.997.107-.775.418-1.305.762-1.604-2.665-.305-5.467-1.334-5.467-5.931 0-1.311.469-2.381 1.236-3.221-.124-.303-.535-1.524.117-3.176 0 0 1.008-.322 3.301 1.23.957-.266 1.983-.399 3.003-.404 1.02.005 2.047.138 3.006.404 2.291-1.552 3.297-1.23 3.297-1.23.653 1.653.242 2.874.118 3.176.77.84 1.235 1.911 1.235 3.221 0 4.609-2.807 5.624-5.479 5.921.43.372.823 1.102.823 2.222v3.293c0 .319.192.694.801.576 4.765-1.589 8.199-6.086 8.199-11.386 0-6.627-5.373-12-12-12z"/>
                </svg>
                {t.hero.viewOnGithub}
              </a>
            </div>
          </div>
          
          <div className="md:w-1/2 flex justify-center">
            <div className="relative w-full max-w-md">
              <div className="absolute -inset-1 bg-gradient-to-r from-blue-500 to-indigo-500 rounded-lg blur opacity-30"></div>
              <div className="relative bg-gray-900 p-2 rounded-lg shadow-xl">
                <img 
                  src="/src/assets/images/pangenome_venn.png" 
                  alt="Pangenome Venn Diagram" 
                  className="w-full h-auto rounded"
                />
                <div className="absolute bottom-0 left-0 right-0 bg-gradient-to-t from-gray-900/90 to-transparent p-4">
                  <p className="text-sm text-gray-300 font-mono">
                    {language === 'zh' ? '泛基因组核心/可变基因组可视化' : 'Pangenome Core/Variable Genome Visualization'}
                  </p>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
}
