import { useI18nStore, translations } from '../lib/i18n';

export function ModulesSection() {
  const { language } = useI18nStore();
  const t = translations[language];

  return (
    <section id="modules" className="py-16 bg-gray-50">
      <div className="container mx-auto px-4">
        <h2 className="text-3xl font-bold text-center mb-12 text-gray-800">
          {t.modules.title}
        </h2>

        {/* 第一部分：基因组分析基础工具 */}
        <div className="mb-16">
          <div className="flex flex-col md:flex-row items-center mb-8">
            <div className="md:w-1/3 mb-6 md:mb-0 md:pr-8">
              <h3 className="text-2xl font-semibold text-blue-800 mb-4">
                {t.modules.part1.title}
              </h3>
              <div className="relative">
                <div className="absolute -inset-1 bg-gradient-to-r from-blue-500 to-indigo-500 rounded-lg blur opacity-30"></div>
                <div className="relative bg-white p-2 rounded-lg shadow-lg">
                  <img 
                    src="/src/assets/images/pangenome_workflow.png" 
                    alt="Pangenome Workflow" 
                    className="w-full h-auto rounded"
                  />
                </div>
              </div>
            </div>
            <div className="md:w-2/3 grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
              <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition-shadow border-l-4 border-blue-500">
                <h4 className="text-lg font-semibold text-gray-800 mb-2">
                  {t.modules.part1.module1.title}
                </h4>
                <p className="text-gray-600">
                  {t.modules.part1.module1.description}
                </p>
              </div>
              <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition-shadow border-l-4 border-blue-500">
                <h4 className="text-lg font-semibold text-gray-800 mb-2">
                  {t.modules.part1.module2.title}
                </h4>
                <p className="text-gray-600">
                  {t.modules.part1.module2.description}
                </p>
              </div>
              <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition-shadow border-l-4 border-blue-500">
                <h4 className="text-lg font-semibold text-gray-800 mb-2">
                  {t.modules.part1.module3.title}
                </h4>
                <p className="text-gray-600">
                  {t.modules.part1.module3.description}
                </p>
              </div>
              <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition-shadow border-l-4 border-blue-500">
                <h4 className="text-lg font-semibold text-gray-800 mb-2">
                  {t.modules.part1.module4.title}
                </h4>
                <p className="text-gray-600">
                  {t.modules.part1.module4.description}
                </p>
              </div>
              <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition-shadow border-l-4 border-blue-500">
                <h4 className="text-lg font-semibold text-gray-800 mb-2">
                  {t.modules.part1.module5.title}
                </h4>
                <p className="text-gray-600">
                  {t.modules.part1.module5.description}
                </p>
              </div>
            </div>
          </div>
        </div>

        {/* 第二部分：群体遗传分析工具 */}
        <div className="mb-16">
          <div className="flex flex-col md:flex-row-reverse items-center mb-8">
            <div className="md:w-1/3 mb-6 md:mb-0 md:pl-8">
              <h3 className="text-2xl font-semibold text-indigo-800 mb-4">
                {t.modules.part2.title}
              </h3>
              <div className="bg-white p-4 rounded-lg shadow-lg border border-indigo-100">
                <div className="aspect-w-16 aspect-h-9">
                  <div className="w-full h-full bg-indigo-50 rounded flex items-center justify-center">
                    <svg className="w-24 h-24 text-indigo-300" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round">
                      <path d="M17 21v-2a4 4 0 0 0-4-4H5a4 4 0 0 0-4 4v2"></path>
                      <circle cx="9" cy="7" r="4"></circle>
                      <path d="M23 21v-2a4 4 0 0 0-3-3.87"></path>
                      <path d="M16 3.13a4 4 0 0 1 0 7.75"></path>
                    </svg>
                  </div>
                </div>
                <div className="mt-4 text-center text-sm text-gray-500 font-mono">
                  {language === 'zh' ? '群体遗传结构分析' : 'Population Genetic Structure Analysis'}
                </div>
              </div>
            </div>
            <div className="md:w-2/3 grid grid-cols-1 md:grid-cols-3 gap-4">
              <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition-shadow border-l-4 border-indigo-500">
                <h4 className="text-lg font-semibold text-gray-800 mb-2">
                  {t.modules.part2.module1.title}
                </h4>
                <p className="text-gray-600">
                  {t.modules.part2.module1.description}
                </p>
              </div>
              <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition-shadow border-l-4 border-indigo-500">
                <h4 className="text-lg font-semibold text-gray-800 mb-2">
                  {t.modules.part2.module2.title}
                </h4>
                <p className="text-gray-600">
                  {t.modules.part2.module2.description}
                </p>
              </div>
              <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition-shadow border-l-4 border-indigo-500">
                <h4 className="text-lg font-semibold text-gray-800 mb-2">
                  {t.modules.part2.module3.title}
                </h4>
                <p className="text-gray-600">
                  {t.modules.part2.module3.description}
                </p>
              </div>
            </div>
          </div>
        </div>

        {/* 第三部分：进化分析工具 */}
        <div className="mb-16">
          <div className="flex flex-col md:flex-row items-center mb-8">
            <div className="md:w-1/3 mb-6 md:mb-0 md:pr-8">
              <h3 className="text-2xl font-semibold text-green-800 mb-4">
                {t.modules.part3.title}
              </h3>
              <div className="relative">
                <div className="absolute -inset-1 bg-gradient-to-r from-green-500 to-teal-500 rounded-lg blur opacity-30"></div>
                <div className="relative bg-white p-2 rounded-lg shadow-lg">
                  <img 
                    src="/src/assets/images/cassava_plant.jpg" 
                    alt="Cassava Plant" 
                    className="w-full h-auto rounded"
                  />
                </div>
              </div>
            </div>
            <div className="md:w-2/3 grid grid-cols-1 md:grid-cols-3 gap-4">
              <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition-shadow border-l-4 border-green-500">
                <h4 className="text-lg font-semibold text-gray-800 mb-2">
                  {t.modules.part3.module1.title}
                </h4>
                <p className="text-gray-600">
                  {t.modules.part3.module1.description}
                </p>
              </div>
              <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition-shadow border-l-4 border-green-500">
                <h4 className="text-lg font-semibold text-gray-800 mb-2">
                  {t.modules.part3.module2.title}
                </h4>
                <p className="text-gray-600">
                  {t.modules.part3.module2.description}
                </p>
              </div>
              <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition-shadow border-l-4 border-green-500">
                <h4 className="text-lg font-semibold text-gray-800 mb-2">
                  {t.modules.part3.module3.title}
                </h4>
                <p className="text-gray-600">
                  {t.modules.part3.module3.description}
                </p>
              </div>
            </div>
          </div>
        </div>

        {/* 第四部分：功能基因组学工具 */}
        <div>
          <div className="flex flex-col md:flex-row-reverse items-center">
            <div className="md:w-1/3 mb-6 md:mb-0 md:pl-8">
              <h3 className="text-2xl font-semibold text-purple-800 mb-4">
                {t.modules.part4.title}
              </h3>
              <div className="bg-white p-4 rounded-lg shadow-lg border border-purple-100">
                <div className="aspect-w-16 aspect-h-9">
                  <div className="w-full h-full bg-purple-50 rounded flex items-center justify-center">
                    <svg className="w-24 h-24 text-purple-300" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round">
                      <polyline points="22 12 18 12 15 21 9 3 6 12 2 12"></polyline>
                    </svg>
                  </div>
                </div>
                <div className="mt-4 text-center text-sm text-gray-500 font-mono">
                  {language === 'zh' ? 'GWAS 分析可视化' : 'GWAS Analysis Visualization'}
                </div>
              </div>
            </div>
            <div className="md:w-2/3 grid grid-cols-1 md:grid-cols-3 gap-4">
              <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition-shadow border-l-4 border-purple-500">
                <h4 className="text-lg font-semibold text-gray-800 mb-2">
                  {t.modules.part4.module1.title}
                </h4>
                <p className="text-gray-600">
                  {t.modules.part4.module1.description}
                </p>
              </div>
              <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition-shadow border-l-4 border-purple-500">
                <h4 className="text-lg font-semibold text-gray-800 mb-2">
                  {t.modules.part4.module2.title}
                </h4>
                <p className="text-gray-600">
                  {t.modules.part4.module2.description}
                </p>
              </div>
              <div className="bg-white p-6 rounded-lg shadow-md hover:shadow-lg transition-shadow border-l-4 border-purple-500">
                <h4 className="text-lg font-semibold text-gray-800 mb-2">
                  {t.modules.part4.module3.title}
                </h4>
                <p className="text-gray-600">
                  {t.modules.part4.module3.description}
                </p>
              </div>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
}
